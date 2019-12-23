'''
$Id: quicklook.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Sun 22 Dec 2019 19:00:22 CET
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

methods for plots of housekeeping and a quicklook at science data
'''
from matplotlib import pyplot as plt
import datetime as dt

def plot_calsource(self,ax=None):
    '''
    plot the calibration source
    '''

    t,v = self.calsource()
    if t is None: return

    dset_shortname = self.dataset_name.split('__')[-1]
    pngname = 'QUBIC_calsource_%s_%s.png' % (dset_shortname,self.obsdate.strftime('%Y%m%d-%H%M%S'))

    tdate = []
    for tstamp in t:
        tdate.append(dt.datetime.fromtimestamp(tstamp))
        
    ttl = 'Calibration Source'

    if ax is None:
        newplot = True
        fontsize = 12
        ttl += '\n'+self.infotext()
        plt.ion()
        fig = plt.figure()
        fig.canvas.set_window_title('plt: Calibration Source for dataset %s' % self.dataset_name)
        fig.suptitle(ttl,fontsize=fontsize)
        ax = fig.add_axes((0.05,0.1,0.9,0.8))
    else:
        newplot = False
        fontsize = 8
        ax.text(0.5,1.0,'Calibration Source',va='bottom',ha='center',fontsize=fontsize,transform=ax.transAxes)
            
    ax.plot(tdate,v)
    ax.set_ylabel('Calibration Source Power / arbitrary units')
    ax.set_xlabel('Date / UT')
    if newplot:
        fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    return ax

def plot_temperatures(self,ax,label,ttl):
    '''
    plot a collection of temperatures
    arguments:
         ax is the plot axis (possibly None)
         label is a dictionary of HK key with sensor label
         ttl is the main title of the plot
    '''
    val = {}
    for sensor in label.keys():
        val[sensor]  = self.get_hk(sensor)

    t = self.get_hk(data='RaspberryDate',hk='EXTERN_HK')
    if t is None: return
    
    tdate = []
    for tstamp in t:
        tdate.append(dt.datetime.fromtimestamp(tstamp))
        
    dset_shortname = self.dataset_name.split('__')[-1]
    pngname = 'QUBIC_%s_%s_%s.png' % (ttl.replace(' ','_'),dset_shortname,self.obsdate.strftime('%Y%m%d-%H%M%S'))

    if ax is None:
        newplot = True
        fontsize = 12
        ttl += '\n'+self.infotext()
        plt.ion()
        fig = plt.figure()
        fig.canvas.set_window_title('plt: %s for dataset %s' % (ttl,self.dataset_name))
        fig.suptitle(ttl,fontsize=fontsize)
        ax = fig.add_axes((0.05,0.1,0.9,0.8))
    else:
        newplot = False
        fontsize = 8
        ax.text(0.5,1.0,ttl,va='bottom',ha='center',fontsize=fontsize,transform=ax.transAxes)
        
        
    for key in val.keys():
        if val[key] is not None:
            plt.plot(tdate,val[key],label=label[key],marker='D')
    ax.set_ylabel('Temperature / K',fontsize=fontsize)
    ax.set_xlabel('Date / UT',fontsize=fontsize)
    ax.legend()
    if newplot:
        fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    return ax
    

def plot_300mKtemperatures(self,ax=None):
    '''
    plot the 300mK temperatures
    '''

    label = {}
    label['AVS47_1_CH2'] = 'TES stage RIRT'
    label['AVS47_1_CH5'] = 'Film breaker'
    label['AVS47_1_CH6'] = '0.3K fridge CH'

    ax = self.plot_temperatures(ax,label,'300mK Temperatures')
    return ax


def plot_1Ktemperatures(self,ax=None):
    '''
    plot the 1K temperatures
    '''
    label = {}
    label['AVS47_1_ch1'] = '1K stage'
    label['AVS47_1_ch3'] = 'M1'
    label['AVS47_1_ch4'] = '1K fridge CH'
    label['AVS47_1_ch7'] = 'M2'
    label['AVS47_2_ch0'] = 'PT2 S2 CH'
    label['AVS47_2_ch2'] = 'Fridge plate MHS'
    label['AVS47_2_ch3'] = '1K stage back'
    label['AVS47_2_ch4'] = '4K shield Cu braids'

    ax = self.plot_temperatures(ax,label,'1K Temperatures')
    return ax
    
    

def quicklook(self,TES=(54,54)):
    '''
    make a page with diagnostic info
    argument: TES is a list of TES to show as examples (one from ASIC1 and one form ASIC2)
    '''

    ttl = 'Diagnostic for %s' % self.dataset_name
    
    ttl += '\n'+self.infotext()
    plt.ion()
    fig = plt.figure(figsize=(10.5,14))
    fig.canvas.set_window_title('plt: %s for dataset %s' % (ttl,self.dataset_name))
    fig.suptitle(ttl)

    # calsource
    ax = fig.add_axes((0.07,0.7,0.4,0.2))
    self.plot_calsource(ax)

    # 300mK temperatures
    ax = fig.add_axes((0.07,0.45,0.4,0.2))
    self.plot_300mKtemperatures(ax)

    # 1K temperatures
    ax = fig.add_axes((0.55,0.45,0.4,0.2))
    self.plot_1Ktemperatures(ax)

    # example timeline from ASIC 1
    ax = fig.add_axes((0.07,0.15,0.4,0.2))
    self.plot_timeline(asic=1,TES=TES[0],ax=ax)

    # example timeline from ASIC 2
    ax = fig.add_axes((0.55,0.15,0.4,0.2))
    self.plot_timeline(asic=2,TES=TES[1],ax=ax)

    pngname = 'QUBIC_quicklook_%s.png' % self.dataset_name
    fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')    
    
    return
