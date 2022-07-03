'''
$Id: qubiciv.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Fri 01 Jul 2022 08:34:34 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

utilities to read/write QUBIC I-V curves
this is a processed data product
'''
import numpy as np
import datetime as dt
import sys,os,time
from astropy.io import fits

'''
assign the comments for the primary header keys
'''
hdr_comment = {}
hdr_comment['TELESCOP'] ='Telescope used for the observation'
hdr_comment['FILETYPE'] = 'File identification'
hdr_comment['DATASET'] = 'QubicStudio dataset name'
hdr_comment['DATE-OBS'] = 'Date data was taken'
hdr_comment['ANALYSIS'] = 'Analysis type used to determine TES quality'
hdr_comment['ANALYSER'] = 'the person who did the analysis'
hdr_comment['ELOG'] = 'link to the elog entry for the data'
hdr_comment['WIKIPAGE'] = 'link to the wiki page where there is more information'
hdr_comment['FILENAME'] = 'name of this file'
hdr_comment['FILEDATE'] = 'UT date this file was created'


    
def assign_header_values(fpobject,elog='',wikipage=''):
    '''
    initialize the primary header
    '''
    filename = '%s_I-V-curves.fits' % fpobject.dataset_name
    hdr = {}
    for key in hdr_comment.keys():
        hdr[key] = None

    datefmt = '%Y-%m-%dT%H:%M:%S'
    hdr['TELESCOP'] = 'QUBIC'
    hdr['FILETYPE'] = 'IV'
    hdr['DATASET'] = fpobject.dataset_name
    hdr['DATE-OBS'] = fpobject.obsdate.strftime(datefmt)
    hdr['ANALYSIS'] = 'I-V curve'
    hdr['ANALYSER'] = 'qubicpack'
    hdr['ELOG'] = elog
    hdr['WIKIPAGE'] = wikipage
    hdr['FILENAME'] = filename
    hdr['FILEDATE'] = dt.datetime.utcnow().strftime(datefmt)
        
    return hdr


def write_iv_fits(fpobject,elog='',wikipage=''):
    '''
    write the IV curves FITS file
    '''
    verbosity = fpobject.verbosity

    prihdr = fits.Header()
    hdrval = assign_header_values(fpobject,elog,wikipage)

    # make a final HDU with additional info because it gets cut off in the primary header
    # this is like a secondary primary header
    secondprimarylist = []
        
    toolong = False
    for key in hdrval.keys():
        val = hdrval[key]
        comment = self.hdr_comment[key]

        fmt = 'A%i' % max((len(val),len(comment)))
        col = fits.Column(name=key, format=fmt, unit='text', array=[val,comment])
        secondprimarylist.append(col)
                
        cardlen = len(val) + len(comment)
        if cardlen>80:
            toolong = True
            comment = ''
        prihdr[key] = (val,comment)
            

    if toolong:
        prihdr['WARNING'] = ('CUTOFF','go to HDU3 for details')
        print('Truncated information can be found in full in the last HDU')

    # make the primary HDU out of the header     
    prihdu = fits.PrimaryHDU(header=prihdr)

    # prepare the hdulist
    hdulist = [prihdu]

    # the I-V curves are in an HDU for each ASIC
    # because the number of sample points is not necessarily the same between ASICs
    nasics = 0
    for asicobj in fpobject.asic_list:
        if asicobj is None: continue
        nasics += 1
        if verbosity>0: print('ASIC %i' % nasics)
        ivcurves = asicobj.best_iv_curve() # tuple with two arrays of shape (128,npts): V,I for each TES

        Vcolumnlist = []
        Icolumnlist = []
        for TESidx in range(128):
            TESstr = '_TES%03i' % (TESidx+1)
            colV = fits.Column(name='V'+TESstr, format='D', unit='ADU', array=ivcurves[0][TESidx,:])
            Vcolumnlist.append(colV)
            colI = fits.Column(name='I'+TESstr, format='D', unit='ADU', array=ivcurves[1][TESidx,:])
            Icolumnlist.append(colI)
                
        cols  = fits.ColDefs(Vcolumnlist)
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.header['ASIC'] = asicobj.asic
        tbhdu.header['IV_AXIS'] = 'V'
        hdulist.append(tbhdu)

        cols  = fits.ColDefs(Icolumnlist)
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.header['ASIC'] = asicobj.asic
        tbhdu.header['IV_AXIS'] = 'I'
        hdulist.append(tbhdu)
  

    # add the final HDU which is like a secondary primary with all the info which may have been cutoff
    cols  = fits.ColDefs(secondprimarylist)
    secprihdu = fits.BinTableHDU.from_columns(cols)
    hdulist.append(secprihdu)

    # write the FITS file
    thdulist = fits.HDUList(hdulist)
    thdulist.writeto(hdrval['FILENAME'],overwrite=True)
    print('I-V curves written to file: %s' % hdrval['FILENAME'])
    return

def read_iv_fits(filename,verbosity=1):
    '''
    read a set of I-V curves from a fits file
    '''
    if not os.path.is_file(filename):
        print('File not found! %s' % filename)
        return None

    hdulist = fits.open(filename)
    nhdu = len(hdulist)

    if nhdu<4:
        print('This does not appear to be an I-V fits file.  Not enough tables.')
        return None

    nasics = (nhdu-2)/2
    nTES = 128*nasics

    hdr = hdulist[0].header
    if 'TELESCOP' not in hdr.keys() or hdr['TELESCOP']!='QUBIC':
        print('This is not a QUBIC fits file')
        return None

    dataset = hdr['DATASET']
    obsdate = hdr['DATE-OBS']
    elog = hdr['ELOG']
    wikipage = hdr['WIKIPAGE']

    infohdu = hdulist[-1]
    for ttype_idx in range(infohdu.header['TFIELDS']):
        key = 'TTYPE%i' % (ttype_idx+1)
        if infohdu.header[key]=='WIKIPAGE':
            wikipage = infohdu.data.field(ttype_idx)[0]
            continue
        if infohdu.header[key]=='ELOG':
            elog = infohdu.data.field(ttype_idx)[0]
            continue
        

    if verbosity>0:
        print('I-V curves from dataset: %s' % dataset)
        print('Measurement date: %s' % obsdate)
        if elog: print('elog link: %s' % elog)
        if wikipage: print('wikipage: %s' % wikipage)

    ivcurves_list = []
    for idx,hdu in enumerate(hdulist):
        if 'ASIC' not in hdu.header.keys(): continue
        asic = hdu.header['ASIC']
        npts = hdu.header['NAXIS2']
        datkey = hdu.header['IV_AXIS']
        if datkey=='V':
            dat = np.zeros(2,128,npts,dtype=np.float64)
            datidx = 0
        else:
            datidx = 1
            
        if verbosity>0:
            print('reading I-V curve:  %i points %s for ASIC %i' % (npts,datkey,asic))
            
        for TESidx in range(128):
            dat[datidx,TESidx,:] = hdu.data.field(TESidx)            
        ivcurves_list.append(dat)

        
            
