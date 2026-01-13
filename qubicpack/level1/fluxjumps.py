'''
$Id: fluxjumps.py
$auth: Belen Costanza <belen@fcaglp.unlp.edu.ar>
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Mon 23 Jan 2023 15:21:44 CET
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

check data for flux jumps, assign flag, or make correction and assign the "corrected" flag

code originally adapted from Jupyter notebook "Automatic_flux_jump_versions" by Belen Costanza
2026-01-12 17:31:20: updated using "jumps_soft.py" in branch bcostanza-dev of qubicsoft
see:  qubic/scripts/Calibration/Flux_jumps/2025-data/jumps_soft.py
'''
import numpy as np
from scipy.signal import argrelextrema, find_peaks, find_peaks_cwt, savgol_filter 
import bottleneck as bn
from sklearn.cluster import DBSCAN

class fluxjumps:

    def __init__(self, thr=[2e5], window_size=300): 

        """ Class for detection of discontinuities in the data 

        Params: 

        thr = threshold or list of thresholds of the haar filter, if a flux it's higher than the treshold, the time sample would be considered as a fux jump candidate
        window_size = size of the bottleneck moving median 

        Return: 

        nc = number of flux jumps
        xc = time sample of the inition of the flux jump
        xcf = time sample of the end of the flux jump

        """

        self.thr = np.array(thr) 
        self.window_size = window_size

    def haar_function(self, todarray):

        tod_haar = np.zeros(todarray.size)
        xf = bn.move_median(todarray, self.window_size)[self.window_size:]
        tod_haar[self.window_size+self.window_size//2:-self.window_size+self.window_size//2] = xf[:-self.window_size] - xf[self.window_size:]

        return tod_haar

    def find_candidates(self, tod_haar, thr):

        number = 0
        jumps = 0
        thr_used = 0

        #iterate over the amplitude thresholds

        for j,thr_val in enumerate(thr):
            if number == 0: #haven't detected any jump yet

                if max(abs(tod_haar)) >= thr_val: #there's a jump                    
                    number += 1
                    thr_used = thr_val
                    jumps = (abs(tod_haar) >= thr_val) #True in the index where there is flux jumps

        return jumps, thr_used

    def clusters(self, todarray, jumps):

        idx = np.arange(len(todarray))
        idx_jumps = idx[jumps]

        if idx_jumps.size > 1:
            clust = DBSCAN(eps=self.window_size//5, min_samples=1).fit(np.reshape(idx_jumps, (len(idx_jumps),1)))
            nc = np.max(clust.labels_)+1
        else:
            nc = 0.
            idx_jumps = 0.
            clust = 0.

        return nc, idx_jumps, clust

    def initial_start_end(self, nc, idx_jumps, tod_haar, thr_used, clust):

        xc = np.zeros(nc, dtype=int)
        xcf = np.zeros(nc, dtype=int)

        for i in range(nc):

            idx_jumps_from_thr = idx_jumps[clust.labels_ == i]
            idx_delta_end_jump = np.where( abs(tod_haar[idx_jumps_from_thr[-1]:]) < thr_used*0.05 )[0][0]
            idx_delta_start_jump = idx_jumps_from_thr[0] - np.where( abs(tod_haar[:idx_jumps_from_thr[0]]) < thr_used*0.05 )[0][-1]
            #idx_delta_start_jump = np.where( tod_haar[:idx_jumps_from_thr[0]] < thr_used*0.05 )[0][-1]
            xc[i] = idx_jumps_from_thr[0] - idx_delta_start_jump
            xcf[i] = idx_jumps_from_thr[-1] + idx_delta_end_jump
            #delta = xcf - xc

        return xc, xcf #delta 

    def unique(self, xc, xcf):

        xc_unique = np.unique(xc)
        xcf_unique = np.unique(xcf)
        nc_unique = len(xc_unique)

        return nc_unique, xc_unique, xcf_unique

    def change_values(self, xc, xcf, max_gap=10):

        xc2 = []
        xcf2 = []

        i = 0
        while i < len(xc):
            xini = xc[i]
            xfin = xcf[i]
            j = i + 1

            # Agrupar mientras estén dentro del margen
            while j < len(xcf) and xc[j] - xfin <= max_gap:
                xfin = xcf[j]
                j += 1

            xc2.append(xini)
            xcf2.append(xfin)
            i = j

        return xc2, xcf2


    def jumps_detection(self, tt, todarray, consec = True, nc_cond=False):

        tod_haar = self.haar_function(todarray) #1. make the haar filter of the raw TOD

        jumps, thr_used= self.find_candidates(tod_haar, self.thr) #2. if the haar filter is higher than a threshold then is a jump (iterate through an array of possible thresholds)

        nc, idx_jumps, clust = self.clusters(todarray, jumps) #3. Cluster the jumps and find the number of jumps detected in every TES

        if nc_cond == True:
            if nc > 11:
                thr = np.array([3e5]) #higher value for the treshold
                tod_haar = self.haar_function(todarray)
                jumps, thr_used = self.find_candidates(tod_haar, thr)
                nc, idx_jumps, clust = self.clusters(todarray, jumps)

        if nc==0:
            xc=0
            xcf=0
            return nc, xc, xcf, thr_used

        xc, xcf = self.initial_start_end(nc, idx_jumps, tod_haar, thr_used, clust) #5. find the beginning and the end of a jump, also the size of the jump
        nc_unique, xc_unique, xcf_unique = self.unique(xc, xcf)

        if consec == True:
            xc_unique, xcf_unique = self.change_values(xc_unique, xcf_unique)
            nc_unique = len(xc_unique)

        return nc_unique, xc_unique, xcf_unique, thr_used
