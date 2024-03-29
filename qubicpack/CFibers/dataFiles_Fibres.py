"""
$Id: AnaFiber.py
$auth: Sophie Henrot-Versille <versille@lal.in2p3.fr>
$created: Mon 14 Aug 2017 1
$updated: Wed 03 Aug 2022 10:10:33 CEST by Steve for compatibility with changes to qubicpack

dictionnary with data (Carbon fibers-like) binary files to be processed
with run_AnaFiber.py and with the corresponding setup of the data taking
period (fll, Voffset for the detectors and AmpFibre for the amplitude
of the signal of the Fiber)

Files on CC-IN2P3 can be found in directory:
  /sps/qubic/Users/Calib/

"""

dataFibre=[
    
#    {"file":"Data_13_07_2017/Pulses 13072017 3/Sums/2017-07-13 143053/sum-asic1-2017.07.13.143053.bin",
#     "asic":"2",
#     "Voffset":"5V",
#     "I_fll":"40",
#     "AmpFibre":"120mV",
#     "minStep":"20000"},
    {"file":"Data_13_07_2017/Pulses 13072017 4/Sums/2017-07-13 143705/sum-asic1-2017.07.13.143705.bin",
     "asic":"2",
     "Voffset":"4.5V",
     "I_fll":"50",
     "AmpFibre":"120mV",
     "minStep":"10000"},
#    {"file":"Data_13_07_2017/Pulses 13072017 5/Sums/2017-07-13 144234/sum-asic1-2017.07.13.144235.bin",
#     "asic":"2",
#     "Voffset":"4.V",
#     "I_fll":"80",
#     "AmpFibre":"120mV",
#     "minStep":"10000"},
    {"file":"Data_13_07_2017/Pulses 13072017 6/Sums/2017-07-13 144658/sum-asic1-2017.07.13.144658.bin",
     "asic":"2",
     "Voffset":"6V",
     "I_fll":"30",
     "AmpFibre":"120mV",
     "minStep":"10000"},
    {"file":"Data_13_07_2017/Pulses 13072017 7/Sums/2017-07-13 144950/sum-asic1-2017.07.13.144950.bin",
     "asic":"2",
     "Voffset":"5.5V",
     "I_fll":"30",
     "AmpFibre":"120mV",
     "minStep":"10000"},
    {"file":"Data_13_07_2017/Pulses 13072017 8/Sums/2017-07-13 152529/sum-asic1-2017.07.13.152529.bin",
     "asic":"1",
     "Voffset":"7V",
     "I_fll":"50",
     "AmpFibre":"120mV",
     "minStep":"10000"},
    {"file":"Data_13_07_2017/Pulses 13072017 9/Sums/2017-07-13 152824/sum-asic1-2017.07.13.152824.bin",
     "asic":"1",
     "Voffset":"6.5V",
     "I_fll":"50",
     "AmpFibre":"120mV",
     "minStep":"10000"},
    {"file":"Data_13_07_2017/Pulses 13072017 10/Sums/2017-07-13 153104/sum-asic1-2017.07.13.153104.bin",
     "asic":"1",
     "Voffset":"6V",
     "I_fll":"50",
     "AmpFibre":"120mV",
     "minStep":"10000"},    
    {"file":"Data_13_07_2017/Pulses 13072017 11/Sums/2017-07-13 153538/sum-asic1-2017.07.13.153538.bin",
     "asic":"1",
     "Voffset":"6V",
     "I_fll":"50",
     "AmpFibre":"20mV",
     "minStep":"10000"},    
    {"file":"Data_13_07_2017/Pulses 13072017 12/Sums/2017-07-13 153821/sum-asic1-2017.07.13.153821.bin",
     "asic":"1",
     "Voffset":"6.5V",
     "I_fll":"50",
     "AmpFibre":"20mV",
     "minStep":"10000"},    
    {"file":"Data_13_07_2017/Pulses 13072017 13/Sums/2017-07-13 154134/sum-asic1-2017.07.13.154134.bin",
     "asic":"1",
     "Voffset":"7V",
     "I_fll":"50",
     "AmpFibre":"20mV",
     "minStep":"10000"},    
    {"file":"Data_13_07_2017/Pulses 13072017 14/Sums/2017-07-13 154434/sum-asic1-2017.07.13.154434.bin",
     "asic":"1",
     "Voffset":"7V",
     "I_fll":"75",
     "AmpFibre":"20mV",
     "minStep":"10000"},    
    {"file":"Data_13_07_2017/Pulses 13072017 15/Sums/2017-07-13 154717/sum-asic1-2017.07.13.154717.bin",
     "asic":"1",
     "Voffset":"6V",
     "I_fll":"25",
     "AmpFibre":"20mV",
     "minStep":"10000"},    
    ]
