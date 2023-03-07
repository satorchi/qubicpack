'''
$Id: raw.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Mon 06 Mar 2023 16:50:10 CET
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

methods for dealing with the RAW data in QubicStudio datasets
'''
import numpy as np

def read_qubicstudio_raw_fits(self,hdu):
    '''read the raw data
    The HDU passed here as the argument should already have been identified as the Raw HDU

    Each TOD is the integration of 100 samples, with some samples
    filtered out (see RAW_MASK).

    The raw data is generated at 2MHz. The data saved in the Raw files
    is a series of 100 samples for each TES, one after another, so we
    get 100 data points at the high data rate for only one TES, and
    then we get the next TES. 
    '''
    self.printmsg('DEBUG: read_qubicstudio_raw_fits object type is %s' % self.__object_type__,verbosity=3)

    extname = hdu.header['EXTNAME'].strip()

    tdata = self.tdata[-1]

    # check which ASIC
    asic = hdu.header['ASIC_NUM']
    if self.asic is None:
        self.asic = asic
    elif self.asic != asic:
        msg='ERROR! ASIC Id does not match: previous=%i, current=%i' % (self.asic, asic)
        tdata['WARNING'].append(msg)
        self.printmsg(msg)
        self.asic = asic
    self.printmsg('Reading raw data for ASIC %i' % asic,verbosity=2)
    
    # initialize hk dictionary
    if extname not in self.hk.keys():
        self.hk[extname] = {}
    self.printmsg("existing keys for hk['%s']: %s" % (extname,self.hk[extname].keys()),verbosity=3)


    # read the number of samples (normally it's 100)
    nsamples = hdu.header['NSAMPLE']
    self.hk[extname]['NSAMPLE'] = nsamples
    
    # read the number of points
    npts = hdu.header['NAXIS2']
    self.hk[extname]['NPTS'] = npts

    # read all the fields
    # CN is the "chronological number" counter tag for each sample so we can see if we lost packets
    # pixelNum is the identifier for the TES whose samples we have for that timestamp
    nfields = hdu.header['TFIELDS']
    for field_idx in range(nfields):
        field_num = field_idx + 1
        key = hdu.header['TTYPE%i' % field_num].strip()
        units = hdu.header['TUNIT%i' % field_num]
        if units=='ms since 1970-01-01T00:00:00': # convert timestamp to seconds instead of milliseconds
            dat = 0.001*hdu.data.field(field_idx)
        else:
            dat = hdu.data.field(field_idx)
        
        if key not in self.hk[extname].keys():
            self.printmsg('storing new %s' % key,verbosity=2)
            self.hk[extname][key] = dat
        else:
            self.printmsg('concatenating %s to pre-existing' % key,verbosity=2)
            self.hk[extname][key] = np.concatenate((self.hk[extname][key],dat))

        
    return True


def exist_raw_data(self):
    '''
    check if we have raw data
    '''
    if not isinstance(self.tdata,list):
        return False
    if len(self.tdata)==0:
        return False

    extname = 'ASIC_RAW'
    if not extname in self.hk.keys():
        return False
    if not 'Raw' in self.hk[extname]:
        return False
    if not isinstance(self.hk[extname]['Raw'],np.ndarray):
        return False
    return True
    

def make_raw_timeline(self,TES=None,axistype='pps'):
    '''
    reconfigure the RAW data array as a timeline
    '''
    if not self.exist_raw_data():
        self.printmsg('No RAW data',verbosity=2)
        return None

    if TES is None:
        self.printmsg('Please enter a valid TES number',verbosity=0)
        return None

    raw = self.hk['ASIC_RAW']['Raw']
    pixnum = self.hk['ASIC_RAW']['pixelNum']
    raw_tstamps = self.timeaxis(datatype='Raw',axistype=axistype)
    
    pixidx = (pixnum==TES)
    nsamples = self.hk['ASIC_RAW']['NSAMPLE']

    ntod = pixidx.sum()
    
    # tstamp_TES = np.zeros((ntod*nsamples))
    # raw_TES = np.zeros((ntod*nsamples))

    tstamp_TES = np.zeros((ntod,nsamples))
    raw_TES = np.zeros((ntod,nsamples))

    # the TOD is the integral of the samples, so the time of the integration is the time of the final sample
    tstamps_samples = -np.flip(np.arange(nsamples)/2e6)
    for idx,tstamp in enumerate(raw_tstamps[pixidx]):
        tstamp_range = tstamp + tstamps_samples
        # istart = idx*nsamples
        # iend = istart + nsamples
        # tstamp_TES[istart:iend] = tstamp_range
        # raw_TES[istart:iend] = raw[pixidx][idx][:]
        tstamp_TES[idx,:] = tstamp_range
        raw_TES[idx,:] = raw[pixidx][idx][:]


    return (tstamp_TES,raw_TES)
