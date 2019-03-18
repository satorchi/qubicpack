"""
$Id: AnaFiber.py
$auth: Sophie Henrot-Versille <versille@lal.in2p3.fr>
$created: Mon 14 Aug 2017 1
set of class to Analyse timeline measurements when Carbon fibers are
pulsed in fron of the TES arrays of QUBIC
"""

from qubicpack import qubicpack as qp
import matplotlib.pyplot as plt
import pandas as pandas
import numpy as numpy


class DataOfOneAsic:
    def __init__(self, asic, timelines=[],tes_blacklist=[],minStep=20000):
        self.asic=asic
        self.timelines = timelines
        self.tes_blacklist=tes_blacklist
        self.minStep=minStep
        self.maxTES=-10
        self.picArray=[]
        self.tesDataObj=[]
        for tes in range(128):
            self.tesDataObj.append(DataOfOneTES(tes,timelines[tes][self.minStep:]))
        self.allSummedData=[]
        self.allMinMax=[]

    def get_ASIC(self):
        return self.asic
    def get_timelines(self):
        return self.timelines
    def get_TESBlackList(self):
        return self.tes_blacklist
    def get_minStep(self):
        return self.minStep
    def get_maxTES(self):
        return self.maxTES
    def get_picArray(self):
        return self.picArray
    def get_allSummedData(self):
        return self.allSummedData
    def get_tesDataObj(self):
        return self.tesDataObj
    def get_allMinMax(self):
        return self.allMinMax

    def determine_TESBlackList(self, doplot=False):
        tes_blacklist=[]
        for tes in range(128):
            if tes == 44:
                plt.plot(self.tesDataObj[tes].get_timeline())
            d,pic_array=self.tesDataObj[tes].compute_summedData(doplot)
            picked=numpy.where(numpy.array(pic_array)<0)
            if numpy.size(picked)<10:
                tes_blacklist.append(tes)
            #if tes == 44:
            #    print numpy.size(picked)
        self.tes_blacklist=tes_blacklist
        return tes_blacklist

    def determine_TESMaxSig(self,doplot=False):
        max_sig=-1000;
        max_tes=-10
        #print 'attention remove first samples a cause de steps..otherwise may need https://github.com/thomasbkahn/step-detect.git'
        for tes in range(128):
            tfromLib2=self.timelines[tes]
            tfromLib=tfromLib2[self.minStep:] 
            if not tes in self.tes_blacklist:
                mean_gliss=pandas.rolling_mean(tfromLib,50)
                mean_gliss=mean_gliss[50:]
                delta_mean=numpy.fabs(mean_gliss.max()-mean_gliss.min())
                if delta_mean>max_sig:
                    max_sig=delta_mean
                    max_tes=tes
                #print delta_mean, max_sig
        print 'chosen tes=',max_tes
        self.maxTES=max_tes
        if doplot:
            data_maxtes=(self.timelines[max_tes])[self.minStep:]
            plt.figure()
            plt.plot(data_maxtes)
            plt.show()

        d,pic_array=self.tesDataObj[max_tes].compute_summedData(doplot)
        self.picArray=pic_array
        return self.maxTES, self.picArray

    def set_TESMaxSig(self,tesMax=1,doplot=False):
        print 'replacing maxTES by ', tesMax
        self.maxTES=tesMax
        if doplot:
            data_maxtes=(self.timelines[max_tes])[self.minStep:]
            plt.figure()
            plt.plot(data_maxtes)
            plt.show()

        d,pic_array=self.tesDataObj[tesMax].compute_summedData(doplot)
        self.picArray=pic_array
        return self.maxTES, self.picArray
    
    def get_allSummed(self,doplot=False):
        allSummedData=[]
        pic_array=self.picArray
        for tes in range(128):
#            if not tes in self.tes_blacklist:
            d=(self.tesDataObj[tes]).sumWithPicArray(pic_array=pic_array)
            allSummedData.append(d)
#            else:
#                allSummedData.append([])
        self.allSummedData=allSummedData
        return allSummedData

    def determine_AllMinMax(self,doplot=False):
        allMinMax=[]
        for tes in range(128):
            d=self.tesDataObj[tes]
            res=d.determine_MinMax()
            allMinMax.append(res)
        self.allMinMax=allMinMax
        return self.allMinMax



    
class DataOfOneTES:
    # data for one TES
    def __init__(self, TES, timeline):
        self.timeline = timeline
        self.TES = TES
        self.summedData=[]
        self.picArray=[]
        self.summedWithPeak=[]
        self.MinMax=[]

    def get_timeline(self):
        return self.timeline
    def get_summedData(self):
        return self.summedData 
    def get_picArray(self):
        return self.pic_Array
    def get_summedWithPeak(self):
        return self.summedWithPeak
    
    def compute_summedData(self,doplot=False):
        datai=self.timeline
        size_data=datai.size
        axe_xdata=numpy.arange(size_data)
        dmean=datai.mean()
        period=0
        max_period=-1000
        min_period=1000
        save_last=0
        pic_array=[]
        count_i=0
        nb_pulse=0
        NbPoints=10
        count_period=0
        summed_signal=[]
        norm_signal=[]
        
        Locked=False
        for i in range(NbPoints,size_data-NbPoints):
            if ( datai[i]<=dmean and datai[i+1]>dmean and datai[i+NbPoints] > datai[i-NbPoints])  :
                count_period+=1
                period = i - save_last
                if period > 200 :
                    if max_period > 0:
                        if Locked==True:
                            nb_pulse+=1
                            mymax=len(summed_signal)
                            for j in range(save_last, mymax):
                                    summed_signal[j-save_last]+=datai[j]
                                    norm_signal[j-save_last]+=1
                        if count_period > 2 and Locked == False:
                            nb_pulse+=1
                            ii=0
                            for j in range(save_last, save_last+period):
                                summed_signal.append(datai[j])
                                norm_signal.append(1)
                                ii+=1
                            Locked=True
                    save_last=i
                    count_i=0
                    if period > max_period: max_period=period
                    if period < min_period: min_period=period
                    pic_array.append(-800000)
            else:
                pic_array.append(1000)
                count_i+=1

        array_signal_Summed = numpy.asarray(summed_signal)
        if doplot:
            plt.figure()
            plt.plot(datai)
            plt.plot(pic_array)
            plt.show()
        array_norm_Summed=numpy.array(norm_signal)
        array_signal_Summed/=array_norm_Summed
        if doplot:
            plt.figure()
            plt.plot(array_signal_Summed-array_signal_Summed.mean(), label=str(self.TES))
            plt.show()
        self.summedData=array_signal_Summed[0:min_period]
        self.picArray=pic_array
        return(array_signal_Summed[0:min_period],pic_array)


    def sumWithPicArray(self,pic_array=[]):
        if numpy.size(pic_array)==0:
            print 'the program will not work ! choose an appropriate pic_array'
        Allsummed=[]
        pic_array_ori=pic_array
        picked=numpy.where(numpy.array(pic_array_ori)<0)
        pic_array=pic_array_ori[picked[0][0]:]
        period_for_sommation=picked[0][1]-picked[0][0]
        summed=[]
        count_period=0
        tfromLib2=self.timeline
        tfromLib=tfromLib2[picked[0][0]:]
        summed=numpy.zeros(period_for_sommation)        
        count_pic=0
        data_pdt_la_periode=numpy.zeros(period_for_sommation)
        for d in range(numpy.size(pic_array)):
            index=d-count_pic
            if index< period_for_sommation and d<picked[0][numpy.size(picked)-1]:
                data_pdt_la_periode[index]=(tfromLib[d])
            if pic_array[d]<0:
                count_pic=d
                if count_pic>0:
                    for i in range(period_for_sommation):
                        summed[i]+=data_pdt_la_periode[i]-data_pdt_la_periode.mean()
                    count_period+=1
                    data_pdt_la_periode=numpy.zeros(period_for_sommation)
        summed/=count_period
        self.summedWithPeak=summed
        return(summed)

    def determine_MinMax(self):
        n,bins,patches=plt.hist(self.summedWithPeak[2:],50)
        self.MinMax=[bins.max(),bins.min()]
        return self.MinMax
        


def plot_Fiber_on_focalPlane(dataToPlot,tagd, xwin=True,figsize=(16,16),color="black",pngname='TES_ARRAY.png'):
    '''
    plot an image of the TES array labeling each pixel
    dataToPlot is an array[asic,pixel]
    '''
    go=qp()
    go.figsize=figsize
    fontsize=figsize[0]
    ttlfontsize=fontsize*1.2
    
    ttl='QUBIC TES array\nASIC1 in blue.  ASIC2 in green.'
    
    nrows=go.pix_grid.shape[0]
    ncols=go.pix_grid.shape[1]
    
    if xwin: plt.ion()
    else: plt.ioff()
    
    count_bad=0
    plt.rc('xtick', labelsize=9)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=9)    # fontsize of the tick labels
    
    fig,ax=plt.subplots(nrows,ncols,figsize=go.figsize)
    if xwin: fig.canvas.set_window_title('plt:  '+ttl)
    fig.suptitle(ttl,fontsize=ttlfontsize)


    TES_translation_table_ASIC1=go.TES2PIX[0]
    TES_translation_table_ASIC2=go.TES2PIX[1]
    
    #WhichTranslat=TES_translation_table_ASIC1
    #if self.asic == "2":
    #        WhichTranslat=TES_translation_table_ASIC2

    for row in range(nrows):
        for col in range(ncols):
            TES=0
            ax[row,col].get_xaxis().set_visible(False)
            ax[row,col].get_yaxis().set_visible(False)
            
            physpix=go.pix_grid[row,col]
            pix_index=physpix-1

            text_y=0.5
            text_x=0.5
            asic_of_pixel=-10
            if physpix==0:
                pix_label=''
                label_colour='black'
                face_colour='black'
            elif physpix in TES_translation_table_ASIC1:
                go.assign_asic(1)
                asic_of_pixel=1
                TES=go.pix2tes(physpix)
                pix_label=str('%i' % TES)
                label_colour='white'
                face_colour='blue'
            elif physpix in TES_translation_table_ASIC2:
                asic_of_pixel=2
                go.assign_asic(2)
                TES=go.pix2tes(physpix)
                pix_label=str('%i' % TES)
                label_colour='white'
                face_colour='green'
            else:
                pix_label='???'
                label_colour='blue'
                face_colour='yellow'
            current=ax[row,col]
            current.set_axis_bgcolor(face_colour)
            current.text(0.5,0.5,pix_label,color="white",fontsize=9, ha='center', va='center', transform=current.transAxes,weight="bold")
            
            if pix_label!='???' and asic_of_pixel>0:
                if tagd[asic_of_pixel-1][int(pix_label)-1] == 1 : # and (int(pix_label)-1 not in self.tes_blacklist):
                    current.set_axis_bgcolor("red")
                    current.text(0.5,0.5,pix_label,color="black",fontsize=9, ha='center', va='center', transform=current.transAxes,weight="bold")
                    count_bad+=1

                alld=dataToPlot[asic_of_pixel-1]
		d=alld[int(pix_label)-1]
	        d=d[2:]
                if numpy.size(d)!=0: 
                    mini=d.min()
                    maxi=d.max()
                    deltai=maxi-mini
                    maxi=maxi+15.*deltai/100.
                    mini=mini-15.*deltai/100.
                    current.set_ylim([mini,maxi])
                    current.plot(d,color=color,linewidth=2.0)

    print 'count_bad',count_bad
    plt.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    if xwin: 
        fig.show()
        plt.show()
    else: plt.close('all')
    return 


 
