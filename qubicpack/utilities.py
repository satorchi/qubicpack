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
import sys,os
import datetime as dt
import matplotlib.pyplot as plt

# on 6 Feb 2018, we reversed the wires for the ASICs
# so now QubicStudio and the dilution fridge use the same ASIC designation
asic_reversal_date = dt.datetime.strptime('2018-02-06 18:00','%Y-%m-%d %H:%M')

# on 27 Feb 2020, qubic-central was finally changed to UTC
qc_utc_date = dt.datetime.strptime('2020-02-27 10:23:33','%Y-%m-%d %H:%M:%S')

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

    # data on qubic-central
    toplevel_dirs.append('/archive')
    
    # data on cc-in2p3
    toplevel_dirs.append('/sps/qubic/Data/Calib-TD')

    # Louise
    toplevel_dirs.append('/home/louisemousset/QUBIC')

    # Steve
    toplevel_dirs.append('/home/work/qubic/data')

    # JCH
    toplevel_dirs.append('/Users/hamilton/Qubic')

    # Martin
    toplevel_dirs.append('/home/martin/QUBIC')

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
        fmt_str = '$%%.%if\\times10^{%%i}$' % nsigfigs
        num_str = fmt_str % (val,expo)
    except:
        num_str = '%.6e' % num

    return num_str
