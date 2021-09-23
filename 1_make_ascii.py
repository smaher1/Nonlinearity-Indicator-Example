#! /usr/bin/env python3.7

# Author: Sean Maher
# Date: April 2020
# This script generates the waveform and spectra files needed to reproduce Figure 3a--3d of Maher et al., (2020)
# The first section windows infrasound data from five stations at Sakurajima Volcano. The MSEED files are available for public download from IRIS. The windowed data are detrended, demeaned and written to Ascii format for plotting.
# The second section performs spectral analysis with a multi-taper cross spectral density function from the Python library NiTime (http://nipy.org/nitime/index.html). Nu_n is computed and the spectra are smoothed using a LOWESS function (Locally Weighted Scatterplot Smoothing) from the module statsmodels (https://www.statsmodels.org/stable/index.html). The smoothed spectra are written to ascii format for plotting with 2_plot_results.py (automatically triggered at the end of this script).

import numpy as np
from obspy import read, UTCDateTime
import scipy.signal
import matplotlib
matplotlib.use("MacOSX")   # fixes plt.show() crashing problem
import matplotlib.pyplot as plt
import glob
import os
import subprocess
import nitime
from statsmodels.nonparametric.smoothers_lowess import lowess


stnlist = ['SVO','KUR','ARI','KOM','HAR']
rangelist = np.array([6263,3408,2400,3487,4525]) # station-vent distance [m]
co = 349
rhonot = 1.225
BETA = 1.201
pref = 20e-6

sfrac=0.0097 # smoothing fraction

# network-coincident STA/LTA trigger-on time and station-specific durations (details in: Matoza, R. S., Fee, D., & Lopez, T. M. (2014). Acoustic characterization of explosion complexity at Sakurajima, Karymsky, and Tungurahua Volcanoes. Seismological Research Letters, 85(6), 1187–1199. https://doi.org/10.1785/0220140110)
trig = '2013-07-22T07:36:05.065'
strstart = trig.replace(':','_').replace('.','_')
durations = np.array([204.3500, 177.9750, 212.9500, 195.3250, 188.4500])

#### first, make .txt ascii files for the waveforms for plotting purposes (fixed start and end times)
for i in range(len(stnlist)):                   # for each station
    msfile = glob.glob('./data/'+str(UTCDateTime(trig).day)+'_07_2013/*'+stnlist[i]+'*')[0]
    st = read(msfile,starttime=UTCDateTime(trig)-5,endtime=UTCDateTime(trig)+35)
    st.detrend('linear')
    st.detrend('demean')
    p = st[0].data
    strstart = trig.replace(':','_')
    strstart = strstart.replace('.','_')
    pfname=stnlist[i]+'_jgr_'+strstart+'.txt'
    print('writing ascii: '+pfname)
    np.savetxt(pfname, p, fmt='%15.10f')


#### second, make .dat ascii files for the power spectra and nu_N
offlist = list()
for i in range(len(stnlist)): # for each station, convert duration to end time
    offlist.append(str(UTCDateTime(trig)+durations[i]))
extrat=20      # extra window time to improve the accuracy of spectral estimates for low frequency components
start = UTCDateTime(trig)-extrat
# Window data to 20 seconds before network-coincident trigger and 20 seconds after station-specific off-time. Detrend and demean the windowed data. Multiply by a Tukey taper that is designed to only affect the amplitudes of the 20 seconds of noise before and after the event window.
for i in range(len(stnlist)):                   # for each station
    end = UTCDateTime(offlist[i])+extrat
    wdur = end-start                     # window duration
    tukshape=2*(extrat/wdur)             # percentage of stream duration that is before/after noise (use for Tukey taper shape)
    file = glob.glob('./data/'+str(start.day)+'_07_2013/*'+stnlist[i]+'*')[0]
    st = read(file,starttime=start,endtime=end)
    st.detrend('linear')
    st.detrend('demean')
    tukey = scipy.signal.tukey(np.size(st[0].data),alpha=tukshape)
    p_tap = st[0].data*tukey
    p_tap2 = p_tap**2
    
    print("     estimating "+stnlist[i]+" spectra (takes a while)")
    s = np.vstack([np.reshape(p_tap,[1,np.size(p_tap)]),np.reshape(p_tap2,[1,np.size(p_tap2)])])
    [f, csd_est] = nitime.algorithms.multi_taper_csd(
                    s,
                    Fs=200,
#                    NW=28*2, #28*2
                    low_bias=False,
                    adaptive=True,
                    sides='onesided')
    psd = np.real(csd_est[0][0])
    quadspec = np.imag(csd_est[1][0])
    
    omg = 2*np.pi*f
    nu_n = -10*np.log10(np.e)*((omg*BETA)/(rhonot*(co**3)))*(quadspec/psd)
    cnu_n = nu_n*rangelist[i]
    
    pfname=stnlist[i]+'_jgr_'+strstart+'.dat'
    print('writing ascii: '+pfname)
    
    # smoothing
    print("     smoothing")
    smooth_f        = lowess(psd, f, is_sorted=True, frac=sfrac, it=0)[:,0]
    smooth_psd      = lowess(psd, f, is_sorted=True, frac=sfrac, it=0)[:,1]
    smooth_quadspec = lowess(quadspec,f, is_sorted=True, frac=sfrac, it=0)[:,1]
    
    smooth_nu_n = -10*np.log10(np.e)*((omg*BETA)/(rhonot*(co**3)))*(smooth_quadspec/smooth_psd)
    smooth_psd_dB   = 10*np.log10(smooth_psd/((20e-6)**2))
    
    smooth_cnu_n = smooth_nu_n*rangelist[i]
    
    outmat = np.hstack([smooth_f.reshape(-1,1), smooth_psd.reshape(-1,1), smooth_psd_dB.reshape(-1,1), smooth_nu_n.reshape(-1,1), smooth_cnu_n.reshape(-1,1)])
    np.savetxt(pfname, outmat, fmt='%15.10f')  # save SVO data

# use this to automatically run the plotting script
subprocess.call(['python','2_plot_results.py'])
