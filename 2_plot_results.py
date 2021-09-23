#! /usr/bin/env python3.7

# Author: Sean Maher
# Date: April 2020

# Reproduces Figures 3a-3d of Maher et al., (2020) using the ascii files created by the script 1_make_ascii.py.


import matplotlib
matplotlib.use("MacOSX")
params = {'text.usetex': True,
    'text.latex.preamble': [r'\usepackage{helvet}', r'\usepackage{amsmath}'],
    'xtick.minor.visible': False, 'xtick.minor.bottom': False, 'grid.linestyle':'--',
    'legend.edgecolor': 'k', 'legend.fancybox': False, 'legend.framealpha': '1'}
matplotlib.rcParams.update(params)
matplotlib.rcParams['lines.markeredgewidth'] = 0.5
matplotlib.rcParams['lines.markeredgecolor'] = 'k'
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal
import glob
import subprocess
import os


############################### functions #################################
def noise_plot_sak():
    noise_low = np.load('5th_avg.npy')
    noise_mid = np.load('50th_avg.npy')
    noise_high_95 = np.load('95th_avg.npy')
    freqs = noise_low[:,0]
    freq2 = noise_mid[:,0]
    freq3 = noise_high_95[:,0]
    psd5 = noise_low[:,1]
    psd50 = noise_mid[:,1]
    psd95 = noise_high_95[:,1]
    return freqs,psd5,psd95,psd50

def freq_bound(psdfile,freqmin=None,freqmax=None):
    freqmin = freqmin or 0
    freqmax = freqmax or 200
    pref=20e-6
    for q in range(1):
        if os.path.getsize(psdfile) > 3:
            [f,Spp,spp2,cospec,qspec,gain,cospec2,phase,Spp_dB,nu_n] = nu_n_calc(psdfile)
            ind_freqmin = (np.abs(f - freqmin)).argmin() # find indices closest to desired max/min freqs
            ind_freqmax = (np.abs(f-freqmax)).argmin()
            f_bounded = f[np.int(ind_freqmin):np.int(ind_freqmax)]
            psd_bounded = Spp[np.int(ind_freqmin):np.int(ind_freqmax)]
            psd_dB_bounded = 10*np.log10(psd_bounded/(pref**2))
            nu_n_bounded = nu_n[np.int(ind_freqmin):np.int(ind_freqmax)]
        else:
            print(plist[q]+' is empty; skipping it')
    return f_bounded,psd_bounded,psd_dB_bounded,nu_n_bounded
##########################################################################




# load waveform files
plist = sorted(glob.glob('*jgr_2013*.txt'))
plist[2], plist[3] = plist[3], plist[2] # switch KOM and KUR so list is ordered by station distance
psdlist = []
for i in range(len(plist)):
    print('waveform files: '+plist[i])
    psdlist.append(plist[i].replace('.txt','.dat'))


# define stations and event time
stnlist = []
for i in range(len(plist)):
    stnlist.append(plist[i][0:3])
trig = '2013-07-22T07:36:05.065'
time = trig.replace(':','_').replace('.','_')
time2 = time.replace('_',':').replace('000Z','')

colorlist = ['r','y','g','b','m']
rangelist = np.array([2400,3408,3487,4525,6263]) # station-vent distance [m]
numlist = np.array([0,2,4,6,8]) # useful for arranging waveforms



# plotting parameters
plt.figure(i+1,figsize=(10,10))
ax1 = plt.subplot(3,1,1)                # all waveforms
ax10 = plt.subplot(3,1,2)               # ARI spectrogram
ax2 = plt.subplot(3,2,5)                # all PSDs
ax3 = plt.subplot(3,2,6)                # all nu_N_tots


default=12
plt.rc('font', size=default)            # controls default text sizes
plt.rc('axes', titlesize=default+2)     # fontsize of the figure title
plt.rc('axes', labelsize=default+2)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=default)      # fontsize of the tick labels
plt.rc('ytick', labelsize=default)      # fontsize of the tick labels
plt.rc('legend', fontsize=default-2)    # legend fontsize
plt.rc('font', family='Helvetica')

for i in range(len(stnlist)):   # for each station
    stn = stnlist[i]
    print("Working on station "+stn)
    p = np.loadtxt(plist[i])
    t = np.arange(np.size(p))/200
    if i == 0:
        pmax = np.max(p)
        [freqs,psd5,psd95,psd50] = noise_plot_sak()
        ax2.fill_between(freqs[1:],psd95[1:],psd5[1:], facecolor='0.8', edgecolor='1', linewidth='0', label='YO noise range',interpolate=True) # plot noise range under PSD
        ax2.plot(freqs, psd50, 'k', ls='--', linewidth=1, label='mid noise approximation')
        ###################### spectrogram ##################
        myNFFT = 256*2
        [Pxx, freqs, bins, im] = ax10.specgram(p, NFFT=myNFFT, Fs=200, noverlap=myNFFT/2, visible=False)
#        [Pxx, freqs, bins, im] = ax10.specgram(p[np.int(15*200):], NFFT=myNFFT, Fs=200, noverlap=myNFFT/2, visible=False)
        ppxx = 10*np.log10(Pxx/((20e-6)**2))
        ax10.pcolormesh(bins, freqs, ppxx, vmin=40, vmax=np.max(ppxx), cmap='jet') # ARI spectrogram
    print("     plotting waveforms")
    ax1.plot(t,(p/pmax)-numlist[i],colorlist[i],linewidth=1.5) # waveforms
    
    ########### plot nitime spectra
    f      = np.loadtxt(psdlist[i],usecols=0)
    psd    = np.loadtxt(psdlist[i],usecols=1)
    psd_dB = np.loadtxt(psdlist[i],usecols=2)
    nu_n   = np.loadtxt(psdlist[i],usecols=3)
    cnu_n  = np.loadtxt(psdlist[i],usecols=4)
    
    print("     plotting spectra")
    ################### PSD #########################################
    ax2.plot(f,psd_dB,colorlist[i],linewidth=2)
    
    ################### nu_N ########################################
    ax3.plot(f,cnu_n,colorlist[i],linewidth=2)

        
ax1.text(t[0]-3, 1.1, r'\textbf{a)}', fontsize=18, fontweight='bold')
ax10.text(t[0]-3, 10, r'\textbf{b)}', fontsize=18, fontweight='bold')
ax2.text(0.14, 125, r'\textbf{c)}', fontsize=18, fontweight='bold')
ax2.text(90, 125, r'\textbf{d)}', fontsize=18, fontweight='bold')

ax1.set_xlim([0,40])
ax1.set_ylim([-9.1,1.1])
ax1.set_xlabel('Time (s)',fontsize=default+2)
ax1.set_ylabel(r'$p/p_{max}$',fontsize=default+2)
ax1.set_title(time2+r';   $p_{max}$ = '+str(np.int(pmax))+' Pa')
ax1.set_yticks([0,-2,-4,-6,-8])
ax1.set_yticklabels((r'ARI',r'HAR',r'KUR',r'KOM',r'SVO'),fontsize=default)
ax1.set_xticks([0,10,20,30,40])
ax1.set_xticklabels((r'0',r'10',r'20',r'30',r'40'),fontsize=default)
    
ax10.set_ylabel('Frequency (Hz)',fontsize=default+2)
ax10.set_xlabel('Time (s)',fontsize=default+2)
ax10.set_xlim([0,40])
ax10.set_ylim([0.1,10])
ax10.set_yticks([0,2,4,6,8,10])
ax10.set_yticklabels((r'0',r'2',r'4',r'6',r'8',r'10'),fontsize=default)
ax10.set_xticks([0,10,20,30,40])
ax10.set_xticklabels((r'0',r'10',r'20',r'30',r'40'),fontsize=default)

ax2.set_xscale('log')
#ax2.set_ylim([0, 120])
ax2.set_xlim([0.1,10])
ax2.set_xticks([0.1,1,10])
ax2.set_yticks([0,20,40,60,80,100,120])
ax2.set_yticklabels([r'0',r'20',r'40',r'60',r'80',r'100',r'120'],fontsize=default)
ax2.set_ylabel(r'PSD (dB re 20 $\mu$Pa$/\sqrt{\text{Hz}}$)',fontsize=default+2)
ax2.set_xlabel('Frequency (Hz)',fontsize=default+2)
ax2.legend(loc='lower left')
ax2.set_xticklabels((r'10$^{-1}$',r'10$^{0}$',r'10$^{1}$'),fontsize=default)


ax3.set_xscale('log')
ax3.set_xlim([0.1,10])
ax3.set_ylim([-0.5,2])
ax3.set_ylabel(r'$\nu_{N_{tot}}$ (dB)')
ax3.set_xlabel('Frequency (Hz)')
ax3.set_yticks([-0.5,0,0.5,1,1.5,2])
ax3.set_yticklabels([-0.5,0,0.5,1,1.5,2])
ax3.set_xticks([0.1,1,10])
ax3.set_xticklabels((r'10$^{-1}$',r'10$^{0}$',r'10$^{1}$'))




plt.subplots_adjust(hspace=0.25,wspace=0.4)
fname='nu_n_result_'+time+'.png'
plt.savefig(fname,dpi=200,bbox_inches='tight')
plt.clf()
subprocess.call(['open',fname])

