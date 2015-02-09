#!/usr/bin/env python

import os,sys
import pyfits
import numpy as np
from scipy import interpolate

import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from scipy.fftpack import fft,ifft
from scipy.interpolate import spline

from functions import gauss_fn

#################################################################
# routine displays nebular emission lines from redshifted galaxy 
# and overlays NIR atmospheric windows
#
# GLW 2015
#################################################################

# usage:
#   NIRsky.py
#      or
#   NIRsky.py galaxy_name 1.9


ph1 = os.getenv("PYTHONHOME1")
if not ph1: raise "PYTHONHOME1 environmental variable not set!"

elines = [["Lyalpha", 1215.67],
          ["CIV",     1550.0],
          ["CIII]",   1909.],
          ["CII]",    2327.],
          ["MgII",    2796.4],
          ["MgII",    2803.5],
          ["NeV",     3326.],
          ["[OII]",   3727.],  # O2
          ["[NeIII]", 3868.7],
          ["Hgamma",  4340.5],  # Hg
          ["[OIII]",  4363.0],  # O31
          ["Hbeta",   4861.3],  # Hb
          ["[OIII]",  4959.0],  # O32
          ["[OIII]",  5007.0],  # O33
          ["[NII]",   6548.1],
          ["Halpha",  6562.8],  # Ha
          ["[NII]",   6583.0],
          ["[SII]",   6717.0],
          ["[SII]",   6731.0],
         ]

# instrument setup
LUCI = 1
HST = 0

# plot variables
y0 = 0
y1 = 1
labels = 1
#lw = 1.0
lw = 2.0
write = 0
sky = 1
trans = 1
interp = 0


#if "--write" in sys.argv:
#   write = 1


if LUCI:
    # K-band
    #w1,w2 = 19500,24500  # Ang
    # K-band zoom
    #w1,w2 = 23000,23600  # Ang
    # NIR
    w1,w2 = 8000,25000  # Ang

elif HST:
    # HST G141
    #w1,w2 = 10500,17900  # Ang
    #w1,w2 = 9500,21500  # Ang
    w1,w2 = 7000,21500  # Ang

wav = np.arange(w1,w2,10.)



obj = None
if len(sys.argv) > 2:
    ################################################
    # print lines available for at specific redshift
    ################################################
    # usage:
    #    NIRsky.py macs1115 1.85
    obj = sys.argv[1] 
    z = float(sys.argv[2])

elif len(sys.argv) > 1:
    ###############################################
    # display atmospheric windows for an instrument
    ###############################################
    # usage:
    #    NIRsky.py __no_object__
    obj = sys.argv[1] 
    z = None

else:
    ############################################
    # print redshift covered by an emission line
    ############################################
    # usage:
    #    NIRsky.py
    eline_ind = 7 # <-- [OII]
    for i,eline in enumerate(elines):
        print i,eline


#fig = plt.figure(figsize=(10,6))
#fig = plt.figure(figsize=(20,12))

#print elines

# Fudge for Wiphu ***REMOVE FOR NORMAL OPERATION***
#i = 0 
#for line,wemit in elines: 
#    print line,wemit
#    if line == "[OII]":    elines.pop(i)
#    elif line == "Hgamma": elines.pop(i)
#    elif line == "[NII]":  elines.pop(i)
#    elif line == "[SII]":  elines.pop(i)
#    i += 1

#print elines

def plot_trans(p):
    transf = "/".join([ph1,"datafiles/sky","mktrans_zm_10_10.dat"])
    L = open(transf).readlines()
    all = [l.split() for l in L if l[0] != "#"]
    wavs = np.array([float(a[0]) for a in all])
    trans = np.array([float(a[1]) for a in all])

    #F = interpolate.interp1d(wavs,trans)
    F = interpolate.InterpolatedUnivariateSpline(wavs,trans)
    Trans = F(wav/1E4)

    #print wavs.shape
    #print wavs
    #p.plot(wavs*1E4,trans,c="0.65",lw=0.1)
    #p.plot(wavs*1E4,trans,c="0.65",lw=0.1)

    if interp:
        p.fill_between(wav,Trans,1,color="0.85")
    else:
        p.fill_between(wavs*1E4,trans,1,color="0.85")

    #p.fill_between(xnew,power_smooth,1,color="0.85")

def plot_sky(p):
    skyf = "/".join([ph1,"datafiles/sky","mk_skybg_zm_10_10_ph.dat"])
    L = open(skyf).readlines()
    all = [l.split() for l in L if l[0] != "#"]
    wavs = np.array([float(a[0]) for a in all])
    #print wavs
    #print wavs[0],wavs[-1]
    #N = 10000.
    #dw = (wavs[-1]-wavs[0])/N
    #nwavs = arange(wavs[0],wavs[-1]+dw,dw)
    #print wavs.shape

    sky = np.array([float(a[1]) for a in all])
    #F = interpolate.interp1d(wavs,sky)
    F = interpolate.InterpolatedUnivariateSpline(wavs,sky)
    Sky = F(wav/10.)



    ax = p.twinx()
    #ax.plot(wav,Sky,c="r",lw=0.1)
    if interp:
        ax.plot(wav,Sky,c="r",lw=0.1)
    else:
        ax.plot(wavs*10,sky,c="r",lw=0.1)
    ax.set_yticks([1,10,100,1000,1E4])
    #ax.set_yticklabels([1,10,100,1000,1E4])
    ax.set_ylim(1,1E4)
    ax.set_yscale("log")
    return ax


def plot_luci_filters(p,text_fontsize=25):
    # http://www.lsw.uni-heidelberg.de/users/agermero/calculator/docu_html/node29.html
    filters = ["z_3002.csv", 
               "J_0403.csv", 
               "H_4302.csv", 
               "K_3902.csv",
               #"Ks_3902.csv",
               "OrderSep_60030.csv",
               "OrderSep_ED763-1.csv",
               ]

    bands = ["z","J","H","K","zJ","HK"]

    for i,filter in enumerate(filters):
        filterf = "/".join([ph1,"datafiles/instruments/LUCI",filter])
        L = open(filterf).readlines()
        all = [l.split() for l in L if l[0] != "#"]
        lines = np.array([float(a[0]) for a in all])
        trans = np.array([float(a[1]) for a in all])

        band = bands[i]
        #if band in ["z","J","H","K"]: c = "g"; h = 0.75; a = 0.85; scale = 1.0
        #else: c = "m"; h = 0.6; a = 0.85; scale = 1.3
        scale = 1.0
        if band in ["z","J","H","K"]: c = "g"; h = 0.75; a = 0.85
        else: c = "m"; h = 0.65; a = 0.85
        p.plot(lines*1E4,trans/scale,alpha=a,lw=2,c=c)
        if labels:
            avg = sum(lines)/len(lines)
            p.text(avg*1E4,h,band,fontsize=text_fontsize,color=c)

def plot_hst_filters(p,text_fontsize=25):
    prefix = "/".join([ph1,"datafiles/instruments/HST/WFC3/"])

    primary = "wfc3_ir_primary_001_syn.fits"
    wfc3 = "wfc3_ir_qe_003_syn.fits"

    grisms  = ["wfc3_ir_g102_src_003_syn.fits","wfc3_ir_g141_src_004_syn.fits"]
    filters = ["wfc3_ir_f105w_004_syn.fits","wfc3_ir_f140w_005_syn.fits"]
    bands   = ["F105W","F140W"]
    disps   = ["G102","G141"]
    
    for i,filter in enumerate(filters):
        grism = grisms[i]

        hdu = pyfits.open(prefix+primary)
        tbdata = hdu[1].data
        # Primary
        wav1 = [t[0] for t in tbdata]
        sen1 = [t[1] for t in tbdata]
        P = interpolate.interp1d(wav1,sen1)
        
        hdu = pyfits.open(prefix+wfc3)
        tbdata = hdu[1].data
        # WFC3
        wav2 = [t[0] for t in tbdata]
        sen2 = [t[1] for t in tbdata]
        W = interpolate.interp1d(wav2,sen2)
        
        hdu = pyfits.open(prefix+grism)
        tbdata = hdu[1].data
        # Grism
        wav3 = [t[0] for t in tbdata]
        sen3 = [t[1] for t in tbdata]

        G = interpolate.interp1d(wav3,sen3)
        
        hdu = pyfits.open(prefix+filter)
        tbdata = hdu[1].data
        # Filter
        wav4 = [t[0] for t in tbdata]
        sen4 = [t[1] for t in tbdata]
        F = interpolate.interp1d(wav4,sen4)
        
        #gwav = np.arange(9701,17906,1)
        gwav = np.arange(int(wav3[0]+1),int(wav3[-1]-1),1)
        
        Primary = P(gwav)
        WFC3 = W(gwav)
        Filter = F(gwav)
        Grism = G(gwav)

        p.plot(gwav,Primary*WFC3*Filter,c="g",lw=lw)
        p.plot(gwav,Primary*WFC3*Grism,c="m",lw=lw)
        #p.plot(wav,Filter)
        
        if labels:
            avg = sum(wav3)/len(wav3)
            p.text(avg-1000,0.65,bands[i],fontsize=text_fontsize,color="g",
                   verticalalignment='center')
            p.text(avg-1000,0.38,disps[i],fontsize=text_fontsize,color="m",
                   verticalalignment='center')




# http://www.gemini.edu/sciops/telescopes-and-sites/observing-condition-constraints/ir-transmission-spectra
# http://www.gemini.edu/sciops/telescopes-and-sites/observing-condition-constraints/ir-background-spectra

#ohf = "OH_skylines.dat"
#ohf = "/".join([ph1,"datafiles/linelists","OH_skylines.dat"])
#L = open(ohf).readlines()
#all = [l.split() for l in L if l[0] != "#"]
#lines = [float(a[0]) for a in all]
#fluxes = [float(a[1]) for a in all]
#for line,flux in zip(lines,fluxes):
#    p.plot([line,line],[0,flux],c="b")

## transparency
#if trans:
#    plot_trans(p)
#
## sky lines
#if sky:
#    plot_sky(p)


#print max(sky)

# convolution
#x = np.arange(len(trans))
#param = 1,0,10,0,0
#gf = gauss_fn(param,0,x)

#A = fft(gf)
#B = fft(trans)
#C = ifft(A*B)
#p.plot(wavs*1E4,C,c="0.65")

fig = plt.figure()
if obj:
    fontsize = 14
    text_fontsize = 20
    mpl.rcParams['xtick.labelsize']=fontsize
    mpl.rcParams['ytick.labelsize']=fontsize

    # transparency
    p = fig.add_subplot(111)
    if trans:
        plot_trans(p)
    # sky lines
    if sky:
        ax = plot_sky(p)
        ax.set_ylabel("Sky Emission",fontsize=fontsize)

    if z:
        for line,wemit in elines:
            wobs = (1+z)*wemit
            p.plot([wobs,wobs],[y0,y1],":",lw=lw,c="b")
            #p.plot([wobs,wobs],[y0,y1],"--",lw=lw,c="b")
            #p.text(wobs+100,0.6,line,fontsize=fontsize,
            if labels:
                p.text(wobs,0.8,line,fontsize=fontsize,
                    rotation='vertical',
                    horizontalalignment='left',
                    verticalalignment='center',)

    if LUCI:
        plot_luci_filters(p)

    elif HST:
        plot_hst_filters(p)

    if labels and z: 
        #p.text(20000,0.95,obj,fontsize=text_fontsize)
        #p.text(20000,0.90,"z = %.3f" % z,fontsize=text_fontsize)    
        p.text(w2-5000,0.95,obj,fontsize=text_fontsize)
        p.text(w2-5000,0.90,"z = %.3f" % z,fontsize=text_fontsize)    
 

    xminorLocator = MultipleLocator(0.1)
    p.yaxis.set_minor_locator(xminorLocator)

    xminorLocator = MultipleLocator(500)
    p.xaxis.set_minor_locator(xminorLocator)

    p.set_xlim(w1,w2)
    p.set_ylim(y0,y1)

    p.set_xlabel("Wavelength ($\AA$)",fontsize=fontsize)
    p.set_ylabel("Sky and Filter Transmission",fontsize=fontsize)

else:
    fontsize=6
    text_fontsize = 8
    mpl.rcParams['xtick.labelsize']=fontsize
    mpl.rcParams['ytick.labelsize']=fontsize

    indices = [0,1,2,3,7,11,13,15]
    plines = [elines[i] for i in indices]
 
    xl = int(np.sqrt(len(plines))+1)

    j = 0
    for i,(line,wemit) in enumerate(plines):
        p = fig.add_subplot(xl,xl,i+1)
        if trans:
            plot_trans(p)
        # sky lines
        if sky:
            ax = plot_sky(p)

        if LUCI:
            obj = "LUCI_filter_curves"
            plot_luci_filters(p,text_fontsize)
            fig.suptitle('LUCI filter curves', fontsize=12)
        
        elif HST:
            obj = "HST_filter_curves"
            plot_hst_filters(p,text_fontsize)
            fig.suptitle('HST filter curves', fontsize=12)

        if labels: 
            #p.text(21000,0.9,line,fontsize=text_fontsize)
            p.text(w2-4000,0.9,line,fontsize=text_fontsize)
        
        xminorLocator = MultipleLocator(0.1)
        p.yaxis.set_minor_locator(xminorLocator)
        xminorLocator = MultipleLocator(1000)
        p.xaxis.set_minor_locator(xminorLocator)
        xmajorLocator = MultipleLocator(5000)
        p.xaxis.set_major_locator(xmajorLocator)
        
        p.set_xlim(w1,w2)
        p.set_ylim(y0,y1)

        ay = p.twiny()
        ay.set_xlim(w1/wemit-1,w2/wemit-1)
        xminorLocator = MultipleLocator(0.1)
        ay.xaxis.set_minor_locator(xminorLocator)

        if j+1 == 1:
            ay.set_xlabel("Redshift (z)",fontsize=fontsize)
        if j+1 == xl or i==xl**2-xl-1:
            p.set_xlabel("Wavelength ($\AA$)",fontsize=fontsize)
        if not i%xl:
            p.set_ylabel("Sky and Filter Transmission",fontsize=fontsize)
        if sky:
            if not (i+1)%xl or i+1 == len(plines):
                ax.set_ylabel("Sky Emission",fontsize=fontsize)

        if not (i+1)%xl:
            j += 1



if write:
    #plt.savefig('%s_NIR.jpg' % obj,orientation='landscape')
    plt.savefig('%s_NIR.pdf' % obj,orientation='landscape')
else:
    plt.show()


