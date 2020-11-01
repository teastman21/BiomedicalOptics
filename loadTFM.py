# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 13:51:29 2020

@author: teast
"""

import numpy as np

def formLoad():
    nDefocusM = np.load("defocusM.npy")
    return nDefocusM

def defocusCalc():
    #global xdf
    xdf = np.empty([M-1],dtype=np.float32)
    fx = (-1 / (2 * dx3) + np.arange(0,M-1,1)*1/L3X)
    #attempt to defocus
    z = 1
    zrange = 1*10**(-3)
    zlist = np.linspace(-zrange,zrange,zs,dtype=np.float32)
    while z < zs:    
        H = np.exp(1j *math.pi * lam1 * zlist[z] * (fx**2))
        H = np.fft.fftshift(H)
        U1=np.fft.fft(np.fft.fftshift(u3x))
        U2 = H * U1
        xdefocus = np.fft.ifftshift(np.fft.ifft(U2))
        xdefocus = np.abs(xdefocus)**2
        xdf = np.vstack((xdf,xdefocus))
        z = z + 1
    return xdf