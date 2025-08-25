'''
$Id: utilities.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Fri 24 May 2019 15:19:42 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

common utilities used in the qubicpack classes
(see also pix2tes.py)
'''
import sys,os,struct
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt

# on 6 Feb 2018, we reversed the wires for the ASICs
# so now QubicStudio and the dilution fridge use the same ASIC designation
asic_reversal_date = dt.datetime.strptime('2018-02-06 18:00','%Y-%m-%d %H:%M')

# on 27 Feb 2020, qubic-central was finally changed to UTC
qc_utc_date = dt.datetime.strptime('2020-02-27 10:23:33','%Y-%m-%d %H:%M:%S')

# on 18 apr 208, we implemented the observation mount used at Alto Chorillos
obsmount_implemented = dt.datetime.strptime('2023-04-18 08:17:10','%Y-%m-%d %H:%M:%S')

# number of pixels in the QUBIC detector matrix per ASIC
NPIXELS = 128


def TES_index(TES):
    '''
    return the index (counting from 0) given the TES number (counting from 1)
    '''
    if TES == 0: # return quietly.  This is used in make_id_focalplane()
        return -1

    if TES>NPIXELS or TES<=0:
        print('TES should have a value between 1 and %i: %i' % (NPIXELS,TES))
        return None
    
    TES_idx=TES-1
    return TES_idx

def ASIC_index(asic):
    '''
    return the asic index (counting from zero) given the asic number
    the asic index is either 0 or 1 which is used for plotting the focal plane

    we do not check for maximum allowed asic_index
    '''
    if asic<1:
        print('asic should have a value greater than 1')
        return None
    asic_idx = asic-1
    remainder = asic_idx % 2
    if remainder == 0:
        return 0
    return 1

def Qubic_DataDir(datadir=None,datafile=None):
    '''
    try to find the user's location for data to be read
    NOTE:  this is different from the qubicpack method assign_datadir which looks for a place to *write* data
    '''
    cwd = os.getcwd() # current working directory

    # make a list of possible directories
    datadirs = []
    toplevel_dirs = []
    
    if datadir is not None:
        toplevel_dirs.append(datadir)

    if 'QUBIC_DATADIR' in os.environ.keys():
        toplevel_dirs.append(os.environ['QUBIC_DATADIR'])

    # Louise
    toplevel_dirs.append('/home/louisemousset/QUBIC')

    # Steve
    toplevel_dirs.append('/home/work/qubic/data')

    # JCH
    toplevel_dirs.append('/Users/hamilton/Qubic')

    # Martin
    toplevel_dirs.append('/home/martin/QUBIC')

    # data on qubic-central
    toplevel_dirs.append('/archive')
    
    # data on cc-in2p3
    toplevel_dirs.append('/sps/qubic/Data/Calib-TD')

    # add the current working directory
    toplevel_dirs.append(cwd)

    for tl_dir in toplevel_dirs:
        for r,f,d in os.walk(tl_dir):
            datadirs.append(r)

    if 'HOME' in os.environ.keys():
        home=os.environ['HOME']
        datadirs.append(home+'/data')

    # now find the first valid data directory
    for datadir in datadirs:
        # if not a valid directory, try the next one
        if not os.path.isdir(datadir): continue

        # if no datafile given, then stop at the first valid directory
        if datafile is None: return datadir
            
        # if a filename is given, check if it's in this directory
        filename = '%s/%s' % (datadir,datafile)
        if os.path.exists(filename): return datadir

    # if we've gone through the whole loop, then we have a problem
    if datafile is None:
        print('ERROR! Could not find a valid data directory')
        return cwd

    print('ERROR! Could not find the directory with that file: %s' % datafile)
    return cwd

def figure_window_title(fig=None,title=None):
    '''
    set a window title for a plot if it's not a Jupyter plot
    '''
    if sys.argv[0].find('ipykernel')>=0: return

    if title is None: ttl = 'plt:'
    else: ttl = 'plt: '+title

    if fig is None: fig = plt.gcf()

    fig.canvas.manager.set_window_title(ttl)
    return
    
def fmt4latex(num,nsigfigs):
    '''
    return a string which formats the number for latex/matplotlib
    '''
    try:
        expo = int(np.log10(num))
        val = num/10**expo
        if val<1:
            expo -= 1
            val *= 10
        if expo==0:
            fmt_str = '$%%.%if$' % nsigfigs
            num_str = fmt_str % val
        else:
            fmt_str = '$%%.%if\\times10^{%%i}$' % nsigfigs
            num_str = fmt_str % (val,expo)
    except:
        num_str = '%.6e' % num

    return num_str

fmt_translation={}
fmt_translation['uint8']   = 'B' 
fmt_translation['int8']    = 'b'
fmt_translation['int16']   = 'h'
fmt_translation['int32']   = 'i'
fmt_translation['int64']   = 'q'
fmt_translation['float32'] = 'f'
fmt_translation['float64'] = 'd'

fmt_reverse_translation = {}
for key in fmt_translation.keys():
    reverse_key = fmt_translation[key]
    fmt_reverse_translation[reverse_key] = key
    

def read_bindat(filename,names=None,fmt=None,STX=0xAA,verbosity=0):
    '''
    read the binary data saved to disk
    '''
    if not os.path.isfile(filename):
        print('ERROR!  File not found: %s' % filename)
        return None

    if names is None:
        print('Please give the names of the data record (comma separated list)')
        return None

    if fmt is None:
        print('Please give the format string (single character per item)')
        return None

    # determine the number of bytes per record entry
    fmt_list = []
    for idx in range(len(fmt)-1):
        fmt_list.append(fmt_reverse_translation[fmt[idx+1]])
    fmt_str = ','.join(fmt_list)    
    rec = np.recarray(names=names,formats=fmt_str,shape=(1))
    nbytes = rec.nbytes

    
    # read the data
    h = open(filename,'rb')
    bindat = h.read()
    h.close()

    # interpret the binary data
    names_list = names.split(',')
    data = {}
    for name in names_list:
        data[name] = []    

    idx = 0
    while idx+nbytes<len(bindat):
        packet = bindat[idx:idx+nbytes]
        dat_list = struct.unpack(fmt,packet)

        if len(dat_list)!=len(names_list):
            print('ERROR:  Incompatible data at byte %i' % idx)
            if verbosity>1: input('enter to continue ')
            idx += 1
            continue

        if dat_list[0]!=STX:
            print('ERROR: Incorrect data at byte %i' % idx)
            if verbosity>1: input('enter to continue ')
            idx += 1
            continue
            

        for datidx,name in enumerate(names_list):
            data[name].append(dat_list[datidx])
            if verbosity>0: print(dat_list)

        idx += nbytes

    for name in data.keys():
        data[name] = np.array(data[name])
        
    return data

def read_obsmount_bindat(filename,verbosity=0):
    '''
    read the binary data acquired from the observation mount and saved to disk
    '''
    rec_fmt = '<Bdd'
    rec_names = 'STX,TIMESTAMP,VALUE'
    return read_bindat(filename,names=rec_names,fmt=rec_fmt,STX=0xAA,verbosity=verbosity)

def read_gps_bindat(filename,verbosity=0):
    '''
    read the binary data acquired from the SimpleRTK GPS system
    '''
    rec_fmt = '<Bdiiiiiiifi'
    rec_names = "STX,timestamp,rpN,rpE,rpD,roll,yaw,pitchIMU,rollIMU,temperature,checksum"
    return read_bindat(filename,names=rec_names,fmt=rec_fmt,STX=0xAA,verbosity=verbosity)

def interpret_rawmask(rawmask_hk,nsamples):
    '''
    interpret the RawMask bit values to sample numbers for the ASIC setting
    '''
    mask = np.zeros(nsamples,dtype=bool)
    for maskidx,bits in enumerate(rawmask_hk):

        if bits==0: continue
        for vector_idx in range(8):
            bitidx = 7 - vector_idx
            bitmask = 2**bitidx
            bitval = (bitmask & bits) >> bitidx
            
            mask[maskidx*8+vector_idx] = bitval
            # print('[%03i] bitmask=%s, bits=%s, bitval=%i' % ((maskidx*8+bitidx),f'{bitmask:08b}',f'{bits:08b}',bitval))
    # nmasked = mask.sum()
    return mask
