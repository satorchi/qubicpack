"""
$Id: ASD.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$auth: Michel Piat
$maintainer: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Wed 05 Jul 2017 17:32:31 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.


plot Amplitude Spectrum Density of time varying TES current
"""
from __future__ import division, print_function
import numpy as np
import sys,os,time
import datetime as dt
import matplotlib.pyplot as plt
from glob import glob
import matplotlib.mlab as mlab
import pickle

from qubicpack.utilities import TES_index
from qubicpack.plot_fp import plot_fp

def plot_ASD(self,TES=None,
             timeline_index=0,
             save=True,
             ax_timeline=None,
             ax_asd=None,
             xwin=True,
             amin=None,amax=None,
             imin=None,imax=None,
             nbins=None,
             indmin=None,indmax=None): #MP
    '''
    plot the Amplitude Spectral Density
    '''

    if not self.exist_timeline_data():
        print('ERROR! No timeline data!')
        return None

    ntimelines=self.ntimelines()
    if timeline_index >= ntimelines:
        print('ERROR! timeline index out of range.  Enter an index between 0 and %i' % (ntimelines-1))
        return None
    
    if TES is None:
        print('Please enter a valid TES number, between 1 and %i' % self.NPIXELS)
        return None
    TES_idx=TES_index(TES)

    if nbins is None:nbins=1
    
    result={}
    result['timeline_index']=timeline_index
    result['TES']=TES
    result['nbins']=nbins
    timeline=self.timeline(TES,timeline_index)
    obsdate=self.tdata[timeline_index]['BEG-OBS']
    result['obsdate']=obsdate

    tinteg=self.tinteg
    if 'INT-TIME' in self.tdata[timeline_index].keys():
        tinteg=self.tdata[timeline_index]['INT-TIME']

    Tbath=self.tdata[timeline_index]['TES_TEMP']
    result['Tbath']=Tbath

    min_bias=self.min_bias
    if 'BIAS_MIN' in self.tdata[timeline_index].keys():
        min_bias=self.tdata[timeline_index]['BIAS_MIN']
    if min_bias is None:min_bias=self.min_bias
    result['min_bias']=min_bias

    max_bias=self.max_bias
    if 'BIAS_MAX' in self.tdata[timeline_index].keys():
        max_bias=self.tdata[timeline_index]['BIAS_MAX']
    if max_bias is None:max_bias=self.max_bias
    result['max_bias']=max_bias

    n_masked=self.n_masked()
    result['n_masked']=n_masked

    current=self.ADU2I(timeline) # uA
    timeline_npts=len(timeline)
    result['timeline_npts']=timeline_npts

    if indmin is None:indmin=0 #MP
    if indmax is None:indmax=timeline_npts-1 #MP

#    bin_npts=timeline_npts//nbins #MP
    bin_npts=(indmax-indmin+1)//nbins #MP
    result['bin_npts']=bin_npts

    if 'NPIXSAMP' in self.tdata[timeline_index].keys():
        npixsampled=self.tdata[timeline_index]['NPIXSAMP']
    else:
        npixsampled=self.NPIXELS
        
    sample_period=self.sample_period()
    time_axis=sample_period*np.arange(timeline_npts)
    
    ttl='Timeline and Amplitude Spectral Density'
    subttl='\nQUBIC Array %s, ASIC %i, TES #%i, T$_\mathrm{bath}$=%.1f mK' % (self.detector_name,self.asic,TES,1000*Tbath)
    pngname='QUBIC_Array-%s_ASIC%i_TES%03i_timeline%03i_AmplitudeSpectralDensity_%s.png'\
             % (self.detector_name,self.asic,TES,timeline_index,obsdate.strftime('%Y%m%dT%H%M%SUTC'))

    pngname_fullpath=self.output_filename(pngname)
    result['pngname']=pngname_fullpath
    
    # setup plot if we haven't done so already
    if xwin: plt.ion()
    else: plt.ioff()
    if ax_timeline is None or ax_asd is None:
        nrows=1
        ncols=2
        fig,axes=plt.subplots(nrows,ncols,sharex=False,sharey=False,figsize=self.figsize)
        ax_timeline=axes[0]
        ax_asd=axes[1]
        if xwin: fig.canvas.set_window_title('plt: '+ttl)
        fig.suptitle(ttl+subttl,fontsize=16)
    result['ax_timeline']=ax_timeline
    result['ax_asd']=ax_asd
    
    txt_x=0.05
    txt_y=0.02

    time_txt=obsdate.strftime('%Y-%m-%m %H:%M:%S UTC')
    time_label='%s Tbath=%.1fmK' % (time_txt,Tbath*1000)
    full_label='%s\nTbath=%.1fmK\nsample period=%.3fmsec\nintegration time=%.1fsec\nnbins=%i\nmin bias=%.2fV\nmax bias=%.2fV\nNpix sampled=%i'\
                % (time_txt,1000*Tbath,1000*sample_period,tinteg,nbins,min_bias,max_bias,npixsampled)
    
    # sampling frequency
    fs = 1.0/self.sample_period()
    
    PSD, freqs = mlab.psd(current[indmin:indmax],
                          Fs = fs,
                          NFFT = bin_npts,
                          window=mlab.window_hanning,
                          detrend='mean')

        
    ax_timeline.cla()
    ax_timeline.plot(time_axis[indmin:indmax],current[indmin:indmax])
    ax_timeline.text(txt_x,txt_y,time_label,transform=ax_timeline.transAxes)
    ax_timeline.set_xlabel('time  /  seconds')
    ax_timeline.set_ylabel('Current  /  $\mu$A')
    if imin is None:imin=min(current)
    if imax is None:imax=max(current)
    ax_timeline.set_ylim((imin,imax))
    if xwin: plt.pause(0.01)

    ASD=np.sqrt(PSD) # in uA
    ax_asd.cla()
    ax_asd.loglog(freqs,ASD)
    boxprops = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax_asd.text(txt_x,txt_y,full_label,transform=ax_asd.transAxes,ha='left',va='bottom',bbox=boxprops)
    ax_asd.set_xlabel('frequency')
    ax_asd.set_ylabel('Amplitude / $\mu$A')
    if amin is None:amin=min(ASD)
    if amax is None:amax=max(ASD)
    ax_asd.set_ylim((amin,amax))
    if xwin: plt.pause(0.01)
        
    if save:
        plt.savefig(pngname_fullpath,format='png',dpi=100,bbox_inches='tight')

    if xwin:plt.show()
    else: plt.close(fig)
    return result

def plot_ASD_all(self,timeline_index=0,imin=None,imax=None,amin=None,amax=None,nbins=None):
    '''
    plot all the ASD for all the TES for a given timeline
    '''
    reslist=[]
    for TES_idx in range(self.NPIXELS):
        TES=TES_idx+1
        result=self.plot_ASD(TES,timeline_index,
                             save=True,ax_timeline=None,ax_asd=None,xwin=False,
                             imin=imin,imax=imax,amin=amin,amax=amax,nbins=nbins)
        reslist.append(result)

    return reslist

def make_ASD_tex_report(self,reslist=None,timeline_index=0):
    '''
    make a tex source file with the report from the ASD measurement
    '''
    if not self.exist_timeline_data():return None

    if reslist is None:
        reslist=self.plot_ASD_all(timeline_index)

    observer=self.observer.replace('<','$<$').replace('>','$>$')
    obsdate=reslist[0]['obsdate']
    Tbath=reslist[0]['Tbath']
    min_bias=reslist[0]['min_bias']
    max_bias=reslist[0]['max_bias']
    timeline_index=reslist[0]['timeline_index']
    timeline_npts=reslist[0]['timeline_npts']
    nbins=reslist[0]['nbins']
    bin_npts=reslist[0]['bin_npts']
    
    texfilename=str('QUBIC_Array-%s_ASIC%i_ASD_%s.tex' % (self.detector_name,self.asic,obsdate.strftime('%Y%m%dT%H%M%SUTC')))
    texfilename_fullpath=self.output_filename(texfilename)
    if texfilename_fullpath is None:
        print('ERROR! Not possible to write tex file.')
        return None
    
    h=open(texfilename_fullpath,'w')
    h.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
    h.write('%%%%% WARNING!  Automatically generated file.  Do not edit! %%%%%\n')
    h.write('%%%%% This file could be overwritten                        %%%%%\n')
    h.write(dt.datetime.utcnow().strftime('%%%%%%%%%% File generated %Y-%m-%d %H:%M:%S UTC                %%%%%%%%%%\n'))
    h.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
    h.write('\\documentclass[a4paper,12pt]{article}\n')
    h.write('\\usepackage{graphicx}\n')
    h.write('\\usepackage{hyperref}\n')
    h.write('\\usepackage{longtable}\n')
    h.write('\\usepackage{setspace}\n')
    h.write('\\newcommand{\\comment}[1]{\n\\begin{minipage}[t]{20ex}\n\\setstretch{0.5}\\flushleft\\noindent\n#1\n\\vspace*{1ex}\n\\end{minipage}}\n')
    
    h.write('\\begin{document}\n')
    h.write('\\begin{center}\n')
    h.write('QUBIC TES Report: Amplitude Spectral Density\\\\\n')
    h.write(obsdate.strftime('data from %Y-%m-%d %H:%M UTC\\\\\n'))
    h.write('compiled by %s\\\\\nusing PyStudio/QubicPack: \\url{https://github.com/satorchi/pystudio}\n' % observer)
    h.write(dt.datetime.utcnow().strftime('this report compiled %Y-%m-%d %H:%M UTC\\\\\n'))
    h.write('\\end{center}\n')

    h.write('\\vspace*{3ex}\n')
    h.write('\\noindent Summary:\n')
    h.write('\\noindent\\begin{itemize}\n')
    h.write('\\item Array %s\n' % self.detector_name)
    h.write('\\item ASIC %i\n' % self.asic)
    if Tbath is None:
        tempstr='unknown'
    else:
        tempstr=str('%.0f mK' % (1000*Tbath))
    h.write('\\item TES physical temperature: %s\n' % tempstr)
    h.write('\\item minimum bias: %.2f V\n' % min_bias)
    h.write('\\item maximum bias: %.2f V\n' % max_bias)
    if min_bias!=max_bias:
        amplitude=0.5*(max_bias-min_bias)
        offset=min_bias+amplitude
        h.write('\\item   equivalent to amplitude %.2f V with offset %.2f V\n' % (amplitude,offset))
    h.write('\\item number of points per timeline: %i\n' % timeline_npts)
    if nbins>1:
        h.write('\\item timeline divided into %i bins with %i points per bin\n' % (nbins,bin_npts))
    h.write('\\end{itemize}\n')
    
    h.write('\n\\vspace*{3ex}\n\\noindent This document includes the following:\n')
    h.write('\\begin{itemize}\n')
    h.write('\\item Plot of all timeline curves in their physical location on the focal plane.\n')
    h.write('\\item Plot of all Spectral Density curves in their physical location on the focal plane.\n')
    h.write('\\item Plot of all the ASD curves.\n')
    h.write('\\end{itemize}\n\\clearpage\n')

    png=str('QUBIC_Array-%s_ASIC%i_timeline_%s.png' % (self.detector_name,self.asic,obsdate.strftime('%Y%m%dT%H%M%SUTC')))
    png_fullpath=self.output_filename(png)
    h.write('\n\\noindent\\includegraphics[width=0.8\\linewidth,clip]{%s}\\\\' % png)

    png=str('QUBIC_Array-%s_ASIC%i_ASD_%s.png' % (self.detector_name,self.asic,obsdate.strftime('%Y%m%dT%H%M%SUTC')))
    png_fullpath=self.output_filename(png)
    h.write('\n\\noindent\\includegraphics[width=0.8\\linewidth,clip]{%s}\n' % png)

    h.write('\\clearpage\n')
    
    for idx,result in enumerate(reslist):
        png=result['pngname']
        h.write('\n\\noindent\\includegraphics[width=0.8\\linewidth,clip]{%s}\\\\' % png)

    
    
    h.write('\n\n\\end{document}\n')
    h.close()
    return texfilename_fullpath
    
def plot_ASD_physical_layout(self,timeline_index=0,xwin=True,amin=None,amax=None,nbins=None):
    '''
    plot the ASD for each TES in it's location in the focal plane
    '''
    if not self.exist_timeline_data():
        print('ERROR! No timeline data!')
        return None

    ntimelines=self.ntimelines()
    if timeline_index >= ntimelines:
        print('ERROR! timeline index out of range.  Enter an index between 0 and %i' % (ntimelines-1))
        return None

    if nbins is None:nbins=1
    
    Tbath=self.tdata[timeline_index]['TES_TEMP']
    obsdate=self.tdata[timeline_index]['BEG-OBS']
    fs = 1.0/self.sample_period()

    pngname=str('QUBIC_Array-%s_ASIC%i_ASD_%s.png' % (self.detector_name,self.asic,obsdate.strftime('%Y%m%dT%H%M%SUTC')))
    pngname_fullpath=self.output_filename(pngname)

    ttl='Amplitude Spectral Density (%s)' % obsdate.strftime('%Y-%m-%d %H:%M')
    subttl='\nQUBIC Array %s, ASIC %i, T$_\mathrm{bath}$=%.1f mK' % (self.detector_name,self.asic,1000*Tbath)

    timeline_npts = self.timeline_array(timeline_index=timeline_index).shape[1]
    psd_npts = timeline_npts // 2
    if np.modf(0.5*timeline_npts)[0] == 0.5:
        psd_npts += 1
    freq_array = np.zeros((self.NPIXELS,psd_npts))
    psd_array = np.zeros((self.NPIXELS,psd_npts))
    for TES_idx in range(self.NPIXELS):

        TES = TES_idx + 1
        timeline=self.timeline(TES,timeline_index)
        current=self.ADU2I(timeline)
        PSD, freqs = mlab.psd(current,
                              Fs = fs,
                              NFFT = timeline_npts//nbins,
                              window=mlab.window_hanning,
                              detrend='mean')
        ASD=np.sqrt(PSD)
        logASD = np.log10(ASD)
        logFreq = np.log10(freqs)

        psd_array[TES_idx,:] = logASD
        freq_array[TES_idx,:] = logFreq

    args = {}
    args['title'] = ttl
    args['subtitle'] = subttl
    args['obsdate'] = obsdate
    args['pngname'] = pngname_fullpath
    key = 'ASIC%i' % self.asic
    args[key] = psd_array
    args['%s x-axis'] = freq_array

    plot_fp(args)
    return args
    
