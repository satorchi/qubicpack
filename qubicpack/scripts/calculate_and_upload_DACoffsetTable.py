#!/usr/bin/env python3
'''
$Id: calculate_and_upload_DACoffsetTable.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Thu 11 Sep 2025 07:38:02 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

calculate the DAC offset table from the given dataset and upload the results back to qubic-central
'''
import sys,os
from qubicpack.qubicfp import qubicfp
from qubichk.utilities import shellcommand

def ADU2offsetDACvalue(ADU,asic_num):
    '''
    convert an ADU value to the offsetDACvalue
    this is slightly modified from the same utility in pystudio where a unsigned int is expected
    '''

    # measured 2026-02-11 and compared to the values found by Michel
    # Michel's values are implemented in pystudio.sequence.get_default_setting

    if asic_num==1:
        coeff = -1.823704185767933e-06
    elif asic_num==2:
        coeff = 1.4861926335974436e-05
    else:
        coeff = 1.4215e-4
        
    offsetDACvalue = coeff*ADU
    return offsetDACvalue


def parseargs(argv):
    '''
    parse some options
    '''
    options = {}
    options['upload'] = True
    options['dataset'] = None
    options['upload directory'] = '.local/share/qubic'
    
    if len(argv)<2:
        print('no arguments')
        return None
    
    for arg in argv[1:]:
        if os.path.isdir(arg):
            options['dataset'] = arg
            continue

        if arg=='--noupload':
            options['upload'] = False
            continue

        if arg.find('upload_dir=')>=0:
            upload_dir = arg.split('=')[1]
            options['upload directory'] = upload_dir
            continue

    return options
        

def help():
    '''
    some help text
    '''
    msg = '\nusage:\n%s <dataset_directory>' % sys.argv[0]
    print(msg)
    return

def upload_offsetTables(file_list,upload_dir=None):
    '''
    upload the offsetDACTables to qubic-central

    NOTE: the server "qubicdl" should be defined in ~/.ssh/config
    see: https://qubic.in2p3.fr/wiki/Operation/ConnectQubic
    '''
    if upload_dir is None:
        upload_dir = '.local/share/qubic'
    cmd = 'scp -p %s qubicdl:%s' % (' '.join(file_list),upload_dir)
    out,err = shellcommand(cmd)
    print(out)
    
    return

def cli():
    '''
    the main program
    '''

    options = parseargs(sys.argv)
    if options is None:
        help()
        return None

    if options['dataset'] is None:
        print('dataset not found or not given')
        return None

    dataset_path = options['dataset']
    a = qubicfp()
    a.read_qubicstudio_dataset(dataset_path)

    # check offset DAC values.  They should all be zero
    file_list = []
    offsetTable_key = 'Offset-DAC-values'
    for idx,asicobj in enumerate(a.asic_list):
        if asicobj is None: continue

        asic_num = idx + 1
        asic_key = 'CONF_ASIC%i' % asic_num
        offsetTable = asicobj.hk[asic_key][offsetTable_key]
        mask = offsetTable!=0
        if mask.sum()!=0:
            print('WARNING!  ASIC %i: Data set has non zero offsets.' % asic_num)

        offsetTable_ADU = asicobj.timeline_array().mean(axis=1)
        offsetTable = ADU2offsetDACvalue(offsetTable_ADU,asic_num)
        ofile_name = 'DAC-Offset-Table_ASIC%02i.txt' % asic_num
        ofile = open(ofile_name,'w')
        for val in offsetTable:
            ofile.write('%15.8e\n' % val)
        ofile.close()
        file_list.append(ofile_name)


    # upload
    if not options['upload']: return
    upload_offsetTables(file_list,options['upload directory'])
    
    return


if __name__=='__main__':
    cli()

    
        
              
