"""
$Id: ana_Fiber_TimeConstant.py
$auth: Sophie Henrot-Versille <versille@lal.in2p3.fr>
$created: Mon 14 Aug 2017 1
$updated: Wed 03 Aug 2022 10:30:51 CEST by Steve for compatibility with qubicpack and python3
re-reading of binary files produced by run_AnaFibre.py which
stores for each file of dataFiles_Fibres and for each TES
the mean template of the signal measured on a pulse of the
Carbon fibers
Will then, for the chosen TES, extract the time constant
(TES+Fibers)
"""

import numpy as numpy
import matplotlib.pyplot as plt
import dataFiles_Fibres as dataFiles
from scipy.optimize import curve_fit  
    
def func(x, a, b, c):
    return a * numpy.exp(-b * x) + c

dataToPlot=numpy.load("plot_dataFibre.bin.npy")


doplot=True
tes_min=0
tes_max=128
listoftes=[8,9,12,15,16,21,25,29,32,34,37,63]
listoffile=[6,7,8,10]
Voff=[]
DeltaV=[]
if doplot:
    get_results=[]
    for i in listoffile:
        plt.figure()
        vec_results=[]
        for tes in listoftes:
            count_fit=0
            dF=dataFiles.dataFibre[i]
            #            if dF["asic"]=="2" and dF["AmpFibre"]=="120mV" and dF["I_fll"]=="50" :                                                                    
            if dF["asic"]=="1" and dF["AmpFibre"]=="20mV" and int(dF["I_fll"])!=75:
                print(i)
                tofit=dataToPlot[i][tes][402:500]
                x=numpy.arange(402,500)
                remove=tofit.min()
                tofit=tofit-remove
                x=x-402
                x=x*0.01
                try:
                    xo,po = curve_fit(func,x,tofit)
                    val=1000./xo[1]
                    res="{0}".format("%d.2" % val)
                    count_fit+=1
                    plt.plot(x,tofit,label="Voffset="+dF["Voffset"]+",I_fll="+dF["I_fll"]+", tau="+str(res)+"ms")
                    plt.xlabel("Time (s)")
                    plt.ylabel("Signal (arbitrary unit)")
                    plt.plot(x,func(x,*xo),'r-')
                    print("|| Voffset=",dF["Voffset"],",I_fll=",dF["I_fll"],'|| ',tes+1,'|| ', res, "||")
                    vec_results.append(res)
                except RuntimeError:
                    print("Error - curve_fit failed")
                    vec_results.append(0)

        print(vec_results)
        get_results.append(vec_results)
        if count_fit==4:
            plt.title("TES="+str(tes+1)+",ASIC="+str(int(dF["asic"])))
            plt.legend()
            #plt.show()
            plt.savefig("compar_sig"+str(tes)+"ASIC1TC.png")
plt.close("all")
plt.figure()
if doplot:
    count=0
    for i in range(numpy.size(listoffile)):
        dF=dataFiles.dataFibre[i]
        plt.plot(listoftes,get_results[i],"o",label="Voffset="+dF["Voffset"]+",I_fll="+dF["I_fll"])
plt.xlabel("tes number ASIC 1")
plt.ylabel("time constant including the C fiber (ms)")
plt.legend(loc=4)
plt.ylim([150,300])
plt.show()

