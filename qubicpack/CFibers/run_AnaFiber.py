"""
$Id: run_AnaFiber.py
$auth: Sophie Henrot-Versille <versille@lal.in2p3.fr>
$created: Mon 14 Aug 2017 1
$updated: Wed 03 Aug 2022 09:48:41 CEST by Steve for compatibility with changes to qubicpack
make use of the AnaFiber set of classes to sum-up the
signal measured for the Carbon fiber pulses and store
them to binary files
"""
import os
import numpy as numpy
import pandas as pandas
import qubicpack.CFibers.AnaFiber as AnaFibre
import matplotlib.pyplot as plt
import qubicpack.CFibers.dataFiles_Fibres as dataFiles
from qubicpack.utilities import Qubic_DataDir
from matplotlib.backends.backend_pdf import PdfPages


doPlot = True

data_dir = Qubic_DataDir(datadir='/sps/qubic/Users/archive/Calib',datafile='Data_13_07_2017')

dataToPlot=[]
allMinMax=[]
allTimeLines=[]
for i in range(numpy.size(dataFiles.dataFibre)):
# pas le 4 ! ? 6 ! ? 
#for i in range(9,10):

    dF=dataFiles.dataFibre[i]
    print(dF["file"])
    filename = os.sep.join([data_dir,dF['file']])
    if not os.path.isfile(filename):
        print('file not found: %s' % filename)
        continue
              
    raw_timelines = AnaFibre.read_bins(filename)
    
    timelines=[]
    timelinessaved=[]
    for tes in range(128):
        test = raw_timelines[tes]
        timelinessaved.append(test)
        pand = pandas.Series(test).rolling(window=1000).mean()
        pand=test[:numpy.size(test)-500]-pand[500:]
        pand=pand[~numpy.isnan(pand)]
        timelines.append(pand[250:])
        
    allTimeLines.append(timelinessaved)

    # create an object DataOfOneAsic
    data=AnaFibre.DataOfOneAsic(dF["asic"], timelines=timelines,minStep=int(dF["minStep"]))
                       
    # then build blackTEstlist:
    bl=data.determine_TESBlackList(doplot=False)
    
    # find the TES with the highest peak-to-peak signal (no glitch)
    # and determine a pic_array 
    a,b=data.determine_TESMaxSig(doplot=False)
    
    if dF["asic"]=="2":
        data.set_TESMaxSig(tesMax=9,doplot=False)
    
    if dF["asic"]=="1":
        data.set_TESMaxSig(tesMax=10,doplot=False)
        
    # get all the summed signal according to the pic_array of
    # the above selected TES 
    allSummed=data.get_allSummed()

    plt.figure()

    if doPlot:
        dataObj=data.get_tesDataObj()
        for tes in data.get_TESBlackList():
            plt.plot(dataObj[tes].get_timeline())
        plt.show()

        plt.figure()
        for tes in range(128):
            plt.plot(allSummed[tes][2:])
        plt.show()
    
    MinMax=data.determine_AllMinMax()
    tMinMax=numpy.transpose(allMinMax)

    if doPlot:
        tes=numpy.arange(128)
        plt.plot(tes,numpy.fabs(tMinMax[1]-tMinMax[0]))
    
    dataToPlot.append(allSummed)
    allMinMax.append(MinMax)


numpy.save("plot_dataFibre.bin",dataToPlot)
numpy.save("minmax_dataFibre.bin",allMinMax)
numpy.save("allTimeLines.bin",allTimeLines)

doplot=False
if doplot:
    for i in range(numpy.size(dataFiles.dataFibre)):
        plt.figure()
        plt.plot(dataToPlot[i][46],label="46")
        plt.plot(dataToPlot[i][45],label="45")
        plt.plot(dataToPlot[i][44],label="44")
        plt.plot(dataToPlot[i][43],label="43")
    #    plt.plot(dataToPlot[i][51],label="51")
        plt.plot(dataToPlot[i][120],label="120")
        plt.plot(dataToPlot[i][124],label="124")
        plt.legend()
        plt.savefig("test_"+str(i)+".png")
        
    for i in range(numpy.size(dataFiles.dataFibre)):
        plt.figure()
        plt.plot(allTimeLines[i][46][100:1000],label="46")
        plt.plot(allTimeLines[i][45][100:1000],label="45")
        plt.plot(allTimeLines[i][44][100:1000],label="44")
        plt.plot(allTimeLines[i][43][100:1000],label="43")
    #    plt.plot(allTimeLines[i][51][100:1000],label="51")
        plt.plot(allTimeLines[i][120][100:1000],label="120")
        plt.plot(allTimeLines[i][124][100:1000],label="124")
        plt.legend()
        plt.savefig("timeline_"+str(i)+".png")

plt.close("all")        
#tes_min=6
#tes_max=7
#tes_min=61
#tes_max=63
#35,38


doplot=False
tes_min=41
tes_max=43
if doplot:
    plt.figure()
    for tes in range(tes_min,tes_max):
        for i in range(numpy.size(dataFiles.dataFibre)):
            dF=dataFiles.dataFibre[i]
            if dF["asic"]=="2" and dF["AmpFibre"]=="120mV" and dF["Voffset"]=="6V" :
                plt.plot(dataToPlot[i][tes]-dataToPlot[i][tes].min(),label="I_fll="+dF["I_fll"])
        plt.title("TES="+str(tes)+",ASIC="+dF["asic"]+",Voffset=6V")
        plt.legend()
        plt.show()
        plt.savefig("compar_sigV6"+str(tes)+".png")
    

    

doplot=False
if doplot:

    plt.figure()
    for tes in range(tes_min,tes_max):
        for i in range(numpy.size(dataFiles.dataFibre)):
            dF=dataFiles.dataFibre[i]
            if dF["asic"]=="2" and dF["AmpFibre"]=="120mV" and dF["I_fll"]=="50" :
                n,bins,patches=plt.hist(dataToPlot[i][tes],50)
                DeltaV.append(numpy.abs(bins.max()-bins.min()))
                Voff.append(float(dF["Voffset"].replace("V","")))
                print(float(dF["Voffset"].replace("V","")), numpy.abs(bins.max()-bins.min()),i)
    plt.figure()
    plt.plot(Voff,DeltaV,label=tes)
    plt.xlim([4,8])
    plt.legend()
    plt.show()

    for tes in range(tes_min,tes_max):
        FindMax=1
        DeltaV=[]
        Voff=[]
        for i in range(numpy.size(dataFiles.dataFibre)):
            dF=dataFiles.dataFibre[i]
            if dF["asic"]=="2" and dF["AmpFibre"]=="120mV" and dF["I_fll"]=="50" :
                if tes==0:
                    print('file',i)
                tMinMax=allMinMax[i][tes]
                Voff.append(float(dF["Voffset"].replace("V","")))
                DeltaV.append(numpy.fabs(tMinMax[0]-tMinMax[1]))
        plt.plot(Voff,DeltaV,label=tes)
    plt.xlim([4,8])
    plt.legend()
    plt.show()
    
