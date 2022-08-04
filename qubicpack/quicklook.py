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
import os
from matplotlib import pyplot as plt
import datetime as dt
import numpy as np
from qubicpack.utilities import qc_utc_date, figure_window_title

def assign_imagename(self,plotname):
    '''
    assign a name to use for the png output
    '''
    dset_shortname = self.dataset_name.split('__')[-1]
    if self.obsdate is None:
        datestr = 'NO-DATE'
    else:
        datestr = self.obsdate.strftime('%Y%m%d-%H%M%S')
        
    imagename = 'QUBIC_%s_%s_%s.png' % (plotname,dset_shortname,datestr)
    return imagename

def plot_calsource(self,ax=None,fontsize=12):
    '''
    plot the calibration source
    '''

    ttl = 'Calibration Source'
    pngname = self.assign_imagename('calsource')
    info = self.calsource_info()
    
    if info is None:
        ttl += ' NO INFORMATION'
    elif info['calsource']['status']!='ON':
        ttl += ' %s' % info['calsource']['status']
    else:
        ttl += ' frequency=%.2fGHz' % info['calsource']['frequency']
    
    if ax is None:
        newplot = True
        ttl = '%s %s' % (self.infotext(),ttl)
        plt.ion()
        fig = plt.figure()
        figure_window_title(fig,'Calibration Source for dataset %s' % self.dataset_name)
        fig.suptitle(ttl,fontsize=fontsize)
        ax = fig.add_axes((0.05,0.1,0.9,0.78))
    else:
        newplot = False
        #ax.text(0.5,1.0,ttl,va='bottom',ha='center',fontsize=fontsize,transform=ax.transAxes)
            
    ax.set_ylabel('Calibration Source Power / arbitrary units',fontsize=fontsize)
    ax.set_xlabel('Date / UT',fontsize=fontsize)
    ax.tick_params(axis='both',labelsize=fontsize)
    ax.text(0.01,1.0,self.calsource_infotext(),va='bottom',ha='left',fontsize=fontsize,transform=ax.transAxes)

    t,v = self.calsource()
    if t is None:
        ax.get_yaxis().set_visible(False)
        ax.get_xaxis().set_visible(False)
        ax.text(0.5,0.5,'No Calsource Data',va='center',ha='center',fontsize=2*fontsize,transform=ax.transAxes)        
    else:
        tdate = []
        for tstamp in t:
            tdate.append(dt.datetime.utcfromtimestamp(tstamp))
        ax.plot(tdate,v)

    if newplot:
        fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    return ax


def plot_temperatures(self,ax,label,ttl,fontsize=12):
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
    
    pngname = self.assign_imagename(ttl.replace(' ','_'))

    if ax is None:
        newplot = True
        ttl += '\n'+self.infotext()
        plt.ion()
        fig = plt.figure()
        figure_window_title(fig,'%s for dataset %s' % (ttl,self.dataset_name))
        fig.suptitle(ttl,fontsize=fontsize)
        ax = fig.add_axes((0.05,0.1,0.9,0.8))
    else:
        newplot = False
        ax.text(0.5,1.0,ttl,va='bottom',ha='center',fontsize=fontsize,transform=ax.transAxes)

    if t is None:
        ax.text(0.5,0.5,'No Temperature Information',va='center',ha='center',fontsize=2*fontsize,transform=ax.transAxes)
        ax.get_yaxis().set_visible(False)
        ax.get_xaxis().set_visible(False)
        if newplot:
            fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
        return ax

    tdate = []
    # qubic-central was changed to UTC on 2020-02-27
    if t[0] > float(qc_utc_date.strftime('%s.%f')):
        for tstamp in t:
            tdate.append(dt.datetime.utcfromtimestamp(tstamp))
    else:
        for tstamp in t:
            tdate.append(dt.datetime.fromtimestamp(tstamp))

    for key in val.keys():
        if val[key] is not None:
            plt.plot(tdate,val[key],label=label[key],marker='D',markersize=0.2*fontsize)
    ax.set_ylabel('Temperature / K',fontsize=fontsize)
    ax.set_xlabel('Date / UT',fontsize=fontsize)
    ax.legend(fontsize=fontsize)
    ax.tick_params(axis='both',labelsize=fontsize)
    if newplot:
        fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    return ax
    

def plot_300mKtemperatures(self,ax=None,fontsize=12):
    '''
    plot the 300mK temperatures
    '''

    label = {}
    label['AVS47_1_CH2'] = 'TES stage RIRT'
    label['AVS47_1_CH5'] = 'Film breaker'
    label['AVS47_1_CH6'] = '0.3K fridge CH'

    ax = self.plot_temperatures(ax,label,'300mK Temperatures',fontsize)
    return ax


def plot_1Ktemperatures(self,ax=None,fontsize=12):
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

    ax = self.plot_temperatures(ax,label,'1K Temperatures',fontsize)
    return ax
    
def plot_switchstatus(self,ax=None,fontsize=12):
    '''
    plot which horn switches are closed
    '''
    v1 = self.get_hk('switch1')
    v2 = self.get_hk('switch2')
    t = self.timeaxis('INTERN_HK')

    ttl = 'Closed Horn Switches'
    pngname = self.assign_imagename(ttl.replace(' ','_'))

    if ax is None:
        newplot = True
        ttl += '\n'+self.infotext()
        plt.ion()
        fig = plt.figure()
        figure_window_title(fig,'%s for dataset %s' % (ttl,self.dataset_name))
        fig.suptitle(ttl,fontsize=fontsize)
        ax = fig.add_axes((0.05,0.1,0.9,0.8))
    else:
        newplot = False
        ax.text(0.5,1.0,ttl,va='bottom',ha='center',fontsize=fontsize,transform=ax.transAxes)


    if t is None or (v1 is None and v2 is None):
        ax.text(0.5,0.5,'No Horn Switch Information',va='center',ha='center',fontsize=2*fontsize,transform=ax.transAxes)
        ax.get_yaxis().set_visible(False)
        ax.get_xaxis().set_visible(False)
    else:
        tdate = []
        for tstamp in t:
            tdate.append(dt.datetime.utcfromtimestamp(tstamp))
        if v1 is not None:
            ax.plot(tdate,v1,marker='D',markersize=0.2*fontsize,ls='none',label='Switch 1 Closed')
        if v2 is not None:
            ax.plot(tdate,v2,marker='D',markersize=0.2*fontsize,ls='none',label='Switch 2 Closed')

        if max(v1)==0 and max(v2)==0:
            ax.text(0.5,0.5,'All horns open',va='center',ha='center',fontsize=2*fontsize,transform=ax.transAxes)
            
        ax.set_ylabel('Horn number',fontsize=fontsize)
        ax.set_xlabel('Date / UT',fontsize=fontsize)
        ax.tick_params(axis='both',labelsize=fontsize)
        ax.set_ylim((-1,65))
        ax.legend(fontsize=fontsize,loc='upper right')
        
    if newplot:
        fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    return ax


def plot_azel(self,ax=None,fontsize=12):
    '''
    plot the azimuth and elevation positions
    '''
    az = self.azimuth()
    el = self.elevation()
    t = self.timeaxis('platform')

    ttl = 'Platform position'
    pngname = self.assign_imagename(ttl.replace(' ','_'))

    if ax is None:
        newplot = True
        ttl += '\n'+self.infotext()
        plt.ion()
        fig = plt.figure()
        figure_window_title(fig,'%s for dataset %s' % (ttl,self.dataset_name))
        fig.suptitle(ttl,fontsize=fontsize)
        ax = fig.add_axes((0.05,0.1,0.9,0.75))
    else:
        newplot = False
        ax.text(0.5,1.0,ttl,va='bottom',ha='center',fontsize=fontsize,transform=ax.transAxes)
        fig = ax.get_figure()

    if t is None:
        ax.text(0.5,0.5,'No Platform information',va='center',ha='center',fontsize=2*fontsize,transform=ax.transAxes)
        ax.get_yaxis().set_visible(False)
        ax.get_xaxis().set_visible(False)
        self.printmsg('No housekeeping information!',verbosity=2)
        if newplot:
            fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
        return ax
        
    
    tdate = []
    for tstamp in t:
        tdate.append(dt.datetime.utcfromtimestamp(tstamp))
    
    if az is None:
        ax.text(0.5,0.5,'No azimuth information',va='center',ha='center',fontsize=2*fontsize,transform=ax.transAxes)
    else:
        ax.plot(tdate,az,marker='D',markersize=0.2*fontsize,ls='none',color='blue',label='Azimuth')
        ax.set_ylabel('Azimuth / Degrees',fontsize=fontsize)    
        ax.legend(loc='lower left',fontsize=fontsize)
        
    if el is None:
        ax.text(0.5,0.4,'No elevation information',va='center',ha='center',fontsize=2*fontsize,transform=ax.transAxes)
    else:
        axel = ax.twinx()
        axel.plot(tdate,el,marker='D',markersize=0.2*fontsize,ls='none',color='red',label='Elevation')
        axel.set_ylabel('Elevation / Degrees',fontsize=fontsize,rotation=270)    
        axel.legend(loc='upper right',fontsize=fontsize)

        
    ax.set_xlabel('Date / UT',fontsize=fontsize)
    ax.tick_params(axis='both',labelsize=fontsize)
    
    if newplot:
        fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    return ax

def plot_hwp(self,ax=None,fontsize=12):
    '''
    plot Half Wave Plate position
    '''

    v = self.get_hk('HWP-Position')
    t = self.timeaxis('hwp')    

    ttl = 'Half Wave Plate position'
    pngname = self.assign_imagename(ttl.replace(' ','_'))

    if ax is None:
        newplot = True
        ttl += '\n'+self.infotext()
        plt.ion()
        fig = plt.figure()
        figure_window_title(fig,'%s for dataset %s' % (ttl,self.dataset_name))
        fig.suptitle(ttl,fontsize=fontsize)
        ax = fig.add_axes((0.05,0.1,0.9,0.75))
    else:
        newplot = False
        ax.text(0.5,1.0,ttl,va='bottom',ha='center',fontsize=fontsize,transform=ax.transAxes)

    if t is None:
        self.printmsg('No housekeeping information!',verbosity=2)
        ax.text(0.5,0.5,'No Half Wave Plate information',va='center',ha='center',fontsize=2*fontsize,transform=ax.transAxes)
        ax.get_yaxis().set_visible(False)
        ax.get_xaxis().set_visible(False)
        if newplot:
            fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
        return ax
    
    tdate = []
    for tstamp in t:
        tdate.append(dt.datetime.utcfromtimestamp(tstamp))

    if v is not None:
        ax.plot(tdate,v,marker='D',markersize=0.2*fontsize,ls='none',label='HWP Position')
        ax.set_ylim(0,8)
        ax.set_ylabel('Position number',fontsize=fontsize)    
        ax.set_xlabel('Date / UT',fontsize=fontsize)
        ax.legend(fontsize=fontsize)
        ax.tick_params(axis='both',labelsize=fontsize)
        
    if v is None or min(v)==255:
        ax.text(0.5,0.5,'No Half Wave Plate information',va='center',ha='center',fontsize=2*fontsize,transform=ax.transAxes)
        ax.get_yaxis().set_visible(False)
        ax.get_xaxis().set_visible(False)

    
    if newplot:
        fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    return ax
    

def quicklook(self,TES=(54,54),
              xwin=True,
              savedir=None,
              filename=None,
              figsize=(10.5,14),
              dpi=100):
    '''
    make a page with diagnostic info
    arguments: 
      - TES       : a list of TES to show as examples (one from ASIC1 and one form ASIC2)
      - xwin      : plot to screen if True, otherwise, just make a file
      - savedir   : directory name for the image file
      - filename  : override the default filename
      - figsize   : override the default figure size.  Note that the intention is for quicklook to fit on A4 size
      - dpi       : override the default dots per inch for the saved image file
    '''
    if not self.exist_data(): return None
    
    ttl = 'Diagnostic for %s' % self.dataset_name
    
    #ttl += '\n'+self.infotext()
    if xwin:
        plt.ion()
    else:
        plt.close('all')
        plt.ioff()
        
    fig = plt.figure(figsize=figsize)
    scalesize = np.sqrt(np.sum(np.array(figsize)**2))
    
    if xwin: figure_window_title(fig,'%s for dataset %s' % (ttl,self.dataset_name))
    fig.suptitle(ttl,fontsize=10)

    fontsize = 0.3*scalesize
    width = 0.37
    height = 0.14
    vspacing = 0.18
    hspacing = 0.51

    hpos1 = 0.07
    hpos2 = hpos1 + hspacing

    vpos = 0.80
    # calsource
    ax = fig.add_axes((hpos1,vpos,width,height))
    self.plot_calsource(ax,fontsize=fontsize)

    # platform position
    ax = fig.add_axes((hpos2,vpos,width,height))
    self.plot_azel(ax,fontsize=fontsize)

    vpos -= vspacing
    # 300mK temperatures
    ax = fig.add_axes((hpos1,vpos,width,height))
    self.plot_300mKtemperatures(ax,fontsize=fontsize)

    # 1K temperatures
    ax = fig.add_axes((hpos2,vpos,width,height))
    self.plot_1Ktemperatures(ax,fontsize=fontsize)

    vpos -= vspacing
    # horn switches activated
    ax = fig.add_axes((hpos1,vpos,width,height))
    self.plot_switchstatus(ax,fontsize=fontsize)

    ax = fig.add_axes((hpos2,vpos,width,height))
    self.plot_hwp(ax,fontsize=fontsize)

    vpos -= vspacing
    # PPS diagnostic
    ax = fig.add_axes((hpos1,vpos,width,height))
    self.plot_pps_nsamples(hk='platform',ax=ax,fontsize=fontsize)

    ax = fig.add_axes((hpos2,vpos,width,height))
    self.plot_pps_nsamples(hk='sci',asic=1,ax=ax,fontsize=fontsize)

    vpos -= vspacing
    # example timeline from ASIC 1
    ax = fig.add_axes((hpos1,vpos,width,height))
    f = self.plot_timeline(asic=1,TES=TES[0],ax=ax,fontsize=fontsize)
    if f is None:
        ax.text(0.5,0.5,'No Science Data for ASIC 1',va='center',ha='center',fontsize=2*fontsize,transform=ax.transAxes)
        ax.get_yaxis().set_visible(False)
        ax.get_xaxis().set_visible(False)


    # example timeline from ASIC 2
    ax = fig.add_axes((hpos2,vpos,width,height))
    f = self.plot_timeline(asic=2,TES=TES[1],ax=ax,fontsize=fontsize)
    if f is None:
        ax.text(0.5,0.5,'No Science Data for ASIC 2',va='center',ha='center',fontsize=2*fontsize,transform=ax.transAxes)
        ax.get_yaxis().set_visible(False)
        ax.get_xaxis().set_visible(False)

    # elapsed time
    delta = self.endobs - self.obsdate
    days = delta.days
    seconds = delta.seconds
    hours = seconds // 3600
    remain_secs = seconds % 3600
    mins = remain_secs // 60
    secs = remain_secs % 60
    usecs = delta.microseconds
    elapsed_time_txt = 'elapsed time:  '
    if days>0:
        elapsed_time_txt += '%i days, ' % days
    elapsed_time_txt += '%i hours, %i minutes, %i.%06i seconds' % (hours,mins,secs,usecs)
    
    vpos -= 0.3*vspacing
    # bottom of page to show elapsed time
    ax = fig.add_axes((hpos1,vpos,width+hspacing,0.2*height))
    ax.text(0.5,0.5,elapsed_time_txt,va='center',ha='center',fontsize=2*fontsize,transform=ax.transAxes)
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)

    if filename is None:
        imagename = 'QUBIC_quicklook_%s.png' % self.dataset_name
    else:
        imagename = filename
        
    if savedir is not None and os.path.isdir(savedir):
        image_fullname = '%s/%s' % (savedir,imagename)
    else:
        image_fullname = imagename
        
    fig.savefig(image_fullname,dpi=dpi,bbox_inches='tight')    

    if not xwin: plt.close(fig)
    
    return image_fullname
