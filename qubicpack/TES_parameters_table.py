#!/usr/bin/env python
'''
$Id: TES_parameters_table.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Fri 15 Dec 2017 14:53:40 CET
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

Make a table which is a list of dictionary values giving the various parameters for each TES detector

Probably the most useful for your scripts is
   TESparameters_dictionary=TEStable_readfile(filename) : read a text version of the TES parameters table

'''
from __future__ import division, print_function
import sys,os,time
import datetime as dt
from qubicpack import qubicpack as qp
import pickle

def read_txtcolumns(filename):
    if not os.path.exists(filename):
        print('ERROR! file not found: %s' % filename)
        return None
    h=open(filename,'r')
    raw=h.read()
    h.close()

    # remove empty lines, and create a 2d list for the rows & columns
    ret=[]
    lines=raw.split('\n')
    idx=0
    for line in lines:
        if line=='':
            del(lines[idx])
            continue
        val=line.split()
        ret.append(val)
        idx+=1
    return ret

def assign_openloop_results(transdic,detector_name,asic):
    '''
    add the Open Loop test results for ASICn
    '''
    filename=str('%s_openlooptest_asic%i.txt' % (detector_name,asic))
    OpenLoopList=read_txtcolumns(filename)
    if not isinstance(OpenLoopList,list):
        print('WARNING! Could not find the Open Loop test results')
        return transdic

    asic_transdic=[entry for entry in transdic if entry['ASIC'] == asic]
    for val in OpenLoopList:
        TES=int(eval(val[0]))
        if val[1]=='0':
            OL='open'
        elif val[1]=='S':
            OL='supra'
        elif val[1]=='1':
            OL='good'        
        else:
            OL=val[1]

        # find the corresponding entry in the list
        TESlist=[entry['TES'] for entry in asic_transdic]
        npts=len(TESlist)
        gotit=False
        idx=0
        while (not gotit) and (idx<npts):
            _TES=TESlist[idx]
            if _TES==TES:
                gotit=True
                entry_idx=idx
            idx+=1
            
        if gotit:
            entry=asic_transdic[entry_idx]
        else:
            transdic.append({})
            entry=transdic[-1]
        entry['OpenLoop']=OL
    
    return transdic

def assign_carbon_fibre_results(transdic,detector_name):
    '''
    assign the identification of bad pixels from the carbon fibre tests by Sophie H-V
    '''
    filename='%s_carbon-fibre_bad-pixel-list.txt' % detector_name
    CarbonFibreList=read_txtcolumns(filename)
    if not isinstance(CarbonFibreList,list):
        print('WARNING! Could not find the Carbon Fibre test results')
        return transdic

    asic_transdic=[]
    asic_transdic.append([entry for entry in transdic if entry['ASIC'] == 1])
    asic_transdic.append([entry for entry in transdic if entry['ASIC'] == 2])
    for val in CarbonFibreList:
        ASIC=int(eval(val[0]))
        asic_index=ASIC-1
        TES=int(eval(val[1]))

        # find the corresponding entry in the list
        TESlist=[entry['TES'] for entry in asic_transdic[asic_index]]
        npts=len(TESlist)
        gotit=False
        idx=0
        while (not gotit) and (idx<npts):
            _TES=TESlist[idx]
            if _TES==TES:
                gotit=True
                entry_idx=idx
            idx+=1
            
        if gotit:
            entry=asic_transdic[asic_index][entry_idx]
        else:
            transdic.append({})
            entry=transdic[-1]
        entry['CarbonFibre']='bad'

    # all other pixels were found to be good by the Carbon Fibre test
    for entry in transdic:
        keys=entry.keys()
        if not 'CarbonFibre' in keys:
            entry['CarbonFibre']='good'
    
    return transdic

def assign_room_temperature_results(transdic,detector_name):
    '''
    assign the Room Temperature Results
    '''
    RTlist=read_txtcolumns('%s_room_temperature.txt' % detector_name)
    if not isinstance(RTlist,list):
        print('WARNING! Could not find the Room Temperature test results')
        return transdic

    go=qp()
    npts_rt=len(RTlist)
    for val in RTlist:
        entry={}
        PIX=int(eval(val[0]))
        entry['PIX']=PIX
        try:
            entry['R300']=eval(val[1])
        except:
            if val[1]=='HS':
                entry['R300']='open'
            else:
                entry['R300']=val[1]

        if PIX in go.TES2PIX[0,:]:
            asic=1
        elif PIX in go.TES2PIX[1,:]:
            asic=2
        else:
            entry['TES']=-1
            entry['ASIC']=-1
            transdic.append(entry)
            continue

        go.assign_asic(asic)
        TES=go.pix2tes(PIX)
        entry['TES']=TES
        entry['ASIC']=asic
        transdic.append(entry)

    return transdic

def make_translation_table(detector_name):
    '''
    create the TES parameters table from various inputs including:
       room temperature resistance (file <array name>_room_temperature.txt)
       carbon fibre evaluation (file carbon-fibre_bad-pixel-list.txt)
       open loop test (file <array name>_openlooptest_asic<n>.txt)

    Note that the input text files should be in the current working directory

    the return is a list of dictionary entries.
    '''    

    transdic=[]

    transdic=assign_room_temperature_results(transdic,detector_name)
    transdic=assign_openloop_results(transdic,detector_name,1)
    transdic=assign_openloop_results(transdic,detector_name,2)

    # assign the results of the Carbon Fibre tests
    transdic=assign_carbon_fibre_results(transdic,detector_name)
    
    # make sure all entries have all keys
    allkeys=['TES','PIX','ASIC','R300','OpenLoop','CarbonFibre']
    nkeys=len(allkeys)
    for entry in transdic:
        for key in allkeys:
            if not key in entry:
                entry[key]='no data'

    # write a pickle file for future reading
    h=open('TES_translation_table.pickle','w')
    pickle.dump(transdic,h)
    h.close()

    return transdic



def TEStable_header(go1,go2):
    '''
    write some explanatory text in the header of the text version of the TES parameters table
    '''
    header ='# QUBIC table of parameters for the TES array\n'
    header+='#\n'
    header+='# all lines beginning with a hash #, such as this one, are comments\n'
    header+='#\n'
    header+='# The first line gives the keywords to be found for each TES\n'
    header+='# Each block of parameters begins with the keyword "INDEX"\n'
    header+='#\n'
    header+='# R300 is the resistance at 300K in units of Ohms\n'
    header+='# CarbonFibre is the evaluation given by the Carbon Fibre tests (good or bad)\n'
    header+='# OpenLoop is the open loop test\n'
    header+='# IV is the evaluation based on the I-V curve (good or bad)\n'
    header+='# K is the coefficient of T in the P_bath vs T_bath plot\n'
    header+='# T0 is the critical temperature in units of K\n'
    header+='# n is the power law index for the P_bath vs T_bath plot\n'
    header+='# NEP is the phonon Noise Equivalent Power in units of W/sqrt(Hz)\n'
    header+='# G is the conductance in units of W/K\n'
    header+='# gamma is the correction coefficient for calculating NEP (see Perbost thesis Eq. 2.72)\n'
    header+='#\n'
    header+='# The following data is for TES array %s\n' % go1.detector_name
    header+='# the IV and NEP evaluation were measured on %s (ASIC %i) and %s (ASIC %i)\n'\
             % (go1.obsdate.strftime('%Y %B %-d'),go1.asic,go2.obsdate.strftime('%Y %B %-d'),go2.asic)
    header+='#\n'
    header+='FILEID=QUBIC table of TES parameters\n'
    return header


def TEStable_entry(idx,nepresults,go):
    '''
    make the entries to the TES parameters table for a given ASIC
    arguments:
      idx : the index of the first entry 
      transdic : the translation table giving the TES results for different tests (room temp, open loop,...)
      nepresults : the table of NEP results (see temperature_analysis.py)
      go  : qubicpack object which has the TES results at the nominal operation temperture (300mK)
    '''
    txt=''
    for TES_index in range(go.NPIXELS):
        TES=TES_index+1
        txt+='\nINDEX=%i' % idx
        entry_both=go.lookup_TEStable(key='TES',value=TES)
        if not entry_both is None:
            entry=[entry for entry in entry_both if entry['ASIC']==go.asic][0]
            for key in entry.keys():
                txt+='\n%s=%s' % (key,str(entry[key]))
            TES=entry['TES']

        if TES<=0:
            txt+='\nIV=unknown'
        else:
            if go.is_good_iv(TES):
                txt+='\nIV=good'
            else:
                txt+='\nIV=bad'

        nepkeys=['K', 'T0', 'n', 'NEP', 'G', 'gamma']
        for nepidx,nepentry in enumerate(nepresults):
            if nepentry['TES']==TES:
                for key in nepkeys:
                    if not nepentry[key] is None:
                        txt+='\n%s=%.8g' % (key,nepentry[key])
                    else:
                        txt+='\n%s=NaN' % key
                key='DET_NAME'
                txt+='\n%s=%s' % (key,nepentry[key])
                break
        idx+=1
    return idx,txt

def TEStable_writefile(nep1,go1,nep2,go2,version=0.0):
    '''
    create the table of TES parameters
    '''
    filename='QUBIC_TES_parameters.v%.1f.txt' % version
    h=open(filename,'w')
    h.write(TEStable_header(go1,go2))
    
    txt='KEYWORDS='
    if not go1.transdic is None:
        for key in go1.transdic[0].keys():
            txt+='%s;' % key
    else:
        print('WARNING! missing results for Room Temperature Resistance, Carbon Fibre, etc')    
    txt+='IV;K;T0;n;NEP;G;gamma;DET_NAME'
    h.write(txt)

    idx=0
    idx,entry=TEStable_entry(idx,nep1,go1)
    h.write(entry)
    h.write('\n')    
    idx,entry=TEStable_entry(idx,nep2,go2)
    h.write(entry)
    h.write('\n')    
    h.close()
    return

def TEStable_readfile(filename):
    '''
    read the text format table of TES parameters
    '''
    if not isinstance(filename,str) or not os.path.exists(filename):
        print('ERROR! File not found: %s' % str(filename))
        return None

    TEStable=[]

    h=open(filename,'r')
    contents=h.read()
    h.close()
    lines=contents.split('\n')

    # ignore all lines starting with '#'
    valid_lines=[]
    for line in lines:
        if not line.find('#')==0 and not line=='':
            valid_lines.append(line)

    # the first line should identify this as a valid file
    line=valid_lines[0]
    keyword,val=line.split('=')
    if not keyword=='FILEID':
        print('ERROR! This does not appear to be a valid QUBIC TES parameters table')
        return None
    if not val=='QUBIC table of TES parameters':
        print('ERROR! This does not appear to be a valid QUBIC TES parameters table')
        return None
    

    # second line should tell us the keywords
    line=valid_lines[1]
    keyword,val=line.split('=')
    if not keyword=='KEYWORDS':
        print('ERROR! Missing the keyword list!')
        return None
    keywords=val.split(';')
    
    msg='Expect to find the following keywords: '
    for keyword in keywords:
        msg+=keyword+' '
    print(msg+'\n')

    val_keywords=['INDEX','R300', 'ASIC', 'PIX', 'TES', 'NEP','G','T0','G','K','n']
    str_keywords=['OpenLoop', 'CarbonFibre', 'IV', 'DET_NAME']
    
    # process the rest of the file
    del(valid_lines[0])
    del(valid_lines[0])
    idx_counter=-1
    entry={}
    for line in valid_lines:
        keyword,str_val=line.split('=')
        if keyword in str_keywords or str_val=='NaN':
            val=str_val
        else:
            val=eval(str_val)
            
        if keyword=='INDEX':
            TEStable.append(entry) # store the previous entry
            entry={}
            idx=val
            idx_counter+=1
            if idx!=idx_counter:
                print('Weirdness: incompatible indeces.  INDEX=%i, idx_counter=%i' % (idx,idx_counter))

        entry[keyword]=val

    # add the last one
    TEStable.append(entry)        
    # clean up:  the first entry is empty because of the way we did the loop
    del(TEStable[0])
    
    
    return TEStable
