"""
$Id: plot_Fibers_FP.py
$auth: Sophie Henrot-Versille <versille@lal.in2p3.fr>
$created: Mon 14 Aug 2017 1
will search for bad pixels according to the data of the TES
recorded with the Carbon fibers and show them on the focal
plane (adapted from Steve Torchinsky's example in scripts/)
"""

import numpy as numpy
import matplotlib.pyplot as plt
import dataFiles_Fibres as dataFiles
import AnaFiber as AnaFibre
from scipy.optimize import curve_fit  
from matplotlib.backends.backend_pdf import PdfPages

    
def func(x, a, b, c):
    return a * numpy.exp(-b * x) + c

dataToPlot=numpy.load("plot_dataFibre.bin.npy")
allTimeLines=numpy.load("allTimeLines.bin.npy")

data=[]

keep_ifile=0
keep_ifile1=0
keep_ifile2=0
ifile=0
for dF in dataFiles.dataFibre:
    if dF["Voffset"] =="6V" and dF["AmpFibre"]=="120mV" and dF["asic"]=="1" :
        data.append(dataToPlot[ifile])
        keep_ifile=ifile
    ifile+=1
ifile1=0
for dF in dataFiles.dataFibre:
    if dF["Voffset"] =="6V" and dF["AmpFibre"]=="120mV" and dF["asic"]=="2" :
        data.append(dataToPlot[ifile1])
        keep_ifile1=ifile1
    ifile1+=1
ifile2=0
for dF in dataFiles.dataFibre:
    if dF["Voffset"] =="6V" and dF["AmpFibre"]=="20mV" and dF["asic"]=="1" and dF["I_fll"]=="50":
        keep_ifile2=ifile2
    ifile2+=1

ifile=keep_ifile
ifile1=keep_ifile1
ifile2=keep_ifile2


from scipy import signal

# ASIC 2
teslist=[1,24,33,61,69,73,88,93,111]

dF=dataFiles.dataFibre[ifile]
dF1=dataFiles.dataFibre[ifile1]
dF2=dataFiles.dataFibre[ifile2]
count_bad=0
with PdfPages('all_timeline.pdf') as pdf:
    tagd=[]
    for tes in teslist:
        plt.figure()
        plt.plot(allTimeLines[ifile1][tes-1][100:10000])
        plt.title("tes="+str(tes)+", asic="+str(dF1["asic"])+"fibre="+str(dF1["AmpFibre"]))
        plt.xlabel("time (sample number)")
        plt.ylabel("signal (arb. units") 
        pdf.savefig()
        plt.close()

# ASIC 1
teslist=[12, 41, 76, 99, 104, 111, 125, 128]

with PdfPages('all_timelineAsic1_less.pdf') as pdf:
    tagd=[]
    for tes in teslist:
        plt.figure()
        plt.plot(allTimeLines[ifile2][tes-1][100:10000])
        plt.title("tes="+str(tes)+", asic="+str(dF2["asic"])+"fibre="+str(dF2["AmpFibre"]))
        plt.xlabel("time (sample number)")
        plt.ylabel("signal (arb. units") 
        pdf.savefig()
        plt.close()

with PdfPages('all_timelineAsic1.pdf') as pdf:
    tagd=[]
    for tes in teslist:
        plt.figure()
        plt.plot(allTimeLines[ifile][tes-1][100:10000])
        plt.title("tes="+str(tes)+", asic="+str(dF["asic"])+"fibre="+str(dF["AmpFibre"]))
        plt.xlabel("time (sample number)")
        plt.ylabel("signal (arb. units") 
        pdf.savefig()
        plt.close()
        
            
count_bad=0
with PdfPages('all.pdf') as pdf:
    tagd=[]
    for asic in range(0,2):
        valtag=[]
        for tes in range(0,128):
            peakind = signal.find_peaks_cwt(data[asic][tes],numpy.arange(1,150))
            valtagdi=0
            if numpy.size(peakind)>4:
#                print 'found peaks !',asic+1, tes+1
                valtagdi=1
            if data[asic][tes][0] == 0. and data[asic][tes][100] == 0. :
                valtagdi=1
#                print 'found flat',asic+1,tes+1

            # compare par reference au 39
            s=1
            if asic==0: 
                s=numpy.sum(numpy.corrcoef(data[0][39],data[asic][tes]))
            if asic==1: 
                s=numpy.sum(numpy.corrcoef(data[1][7],data[asic][tes]))
            if s < 3.582:
                valtagdi=1
                #print 'signal different de ce qu il devrait etre'
                
            if valtagdi != 10:
                plt.figure()
                plt.plot(data[asic][tes])
                plt.title("tes="+str(tes+1)+", asic="+str(asic+1))
                plt.xlabel("time (sample number)")
                plt.ylabel("signal (arb. units") 
                pdf.savefig()
                plt.close()
                count_bad+=1
                print 'Bad Pixel',asic+1,tes+1
            
            valtag.append(valtagdi)
        tagd.append(valtag)
print "bad from real decompte", count_bad
AnaFibre.plot_Fiber_on_focalPlane(data,tagd,xwin=True,figsize=(16,16),color="yellow",pngname="fibreData1307_6.5_50_120.png")                                                                                  
