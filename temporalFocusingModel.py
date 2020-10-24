# -*- coding: utf-8 -*-
"""
Temporal Focusing Model

Created on Wed Sep 30 21:51:19 2020


@author: Tommy Eastman
"""

# import necessary packages 
import sys
sys.modules[__name__].__dict__.clear()
import numpy as np
import math as math
import matplotlib.pyplot as plt
from multiprocessing import Pool

global x1c
global x2c
global x3c
M = 2**14

def initialize():
    """initialize the np arrays for the coordinate build up
        only run once at beginning!
        """
    global x1c
    global x2c
    global x3c
    global xdf

    x1c = np.array([])
    x2c = np.array([])
    x3c = np.array([])
    xdf = np.empty([M-1])


def incident(lam,thetaDcenter):
    #determine diffraction angle
    m = 1 #diffraction order
    g = 120000 #grating density
    lamC = 800*10**(-9) #center wavelength
    thetaI = np.arcsin(m * lamC * g - np.sin(thetaDcenter))
    thetaD = np.arcsin(m * lam * g - np.sin(thetaI))
    return thetaD

def prop(lamm,thetaD): 
    # Define necessary variables
    f1 = 0.5        #collimating lens focal length
    f2 = 9*10**-3   #objective lens focal length
    lamC = 800*10**(-9) #center wavelength
    g = 120000      #grating density
    thetaDcenter = 0  #diffracted center wavelength 
    m = 1           #diffraction order
    w = 0.002       #beam width
    L = 1.05
    
    global zs
    zs = 30        #defocus slices
    global xdf
    x1c = np.array([])
    x2c = np.array([])
    x3c = np.array([])
    
    #multiplied these two by a hundred due to CPU memory constraints @ output
    dx = L / M
    dy = L / M
    x = np.arange(-L/2, L / 2 - dx + .000000000000000001, dx)
    y = np.arange(-L/2, L / 2 - dy + .000000000000000001, dy)
    k = 2 * math.pi / lamm
    #lam1 = 820 * 10 ** (-9)
    lam1 = lamm
    
    #initialize beam and apply tilt
    u1x = np.exp(-x ** 2 / (2 * w ** 2))* np.exp(-1j * k * (x * np.sin(thetaD)))
    
    u1y = np.exp(-y ** 2 / (2 * w ** 2))
    

    #plot beam after tilt applied
    plotU1X = np.abs(u1x)**2

    #x1c = np.vstack((x1c,plotU1X))
    x1c = plotU1X
    
    #propagate to the focal plane (x)
    L2X = lam1 * f1 / dx
    dx2 = lam1 * f1 / L
    
    
    #y
    L2Y = lam1 * f1 / dx
    dy2 = lam1 * f1 / L
    
    #create coordinates at focal plane
    xfoc1 = -L2X / 2 + np.arange(0, M - 1) * dx2
    yfoc1  = -L2Y / 2 + np.arange(0, M - 1) * dy2
    
    #propagate beam with the fourier transform
    u2x = 1/(1j * lam1 * f1)*(np.fft.ifftshift(np.fft.fft(np.fft.fftshift(u1x))))*dx
    u2y = 1/(1j * lam1 * f1)*(np.fft.ifftshift(np.fft.fft(np.fft.fftshift(u1y))))*dy
    #plot beam at fourier plane
  
    plotU2X = np.abs(u2x)**2
    #x2c = np.vstack((x2c,plotU2X))
    x2c = plotU2X
        
    #propagate to output focal plane
    L3X = lam1 * f2 / dx2
    dx3 = lam1 * f2 / L2X
    
    L3Y = lam1 * f2 / dy2
    dy3 = lam1 * f2 / L2Y
     
    #create coordinates at output focal plane
    xfoc2 = -L3X / 2 + np.arange(0, M - 1) * dx3
    yfoc2 = -L3Y / 2 + np.arange(0, M - 1) * dy3
    
    #propagate beam with fourier transform
    u3x = 1/(1j * lam1 * f2)*(np.fft.ifftshift(np.fft.fft(np.fft.fftshift(u2x))))*dx2 
    u3y = 1/(1j * lam1 * f2)*(np.fft.ifftshift(np.fft.fft(np.fft.fftshift(u2y))))*dy2 
    #plot beam at output plane
    
    plotU3X = np.abs(u3x)**2

    #x3c = np.vstack((x3c,plotU3X))
    x3c = plotU3X
    
    fx = (-1 / (2 * dx3) + np.arange(0,M-1,1)*1/L3X)
    #attempt to defocus
    z = 1
    zrange = 1*10**(-3)
    zlist = np.linspace(-zrange,zrange,zs)
    while z < zs:    
        H = np.exp(1j *math.pi * lam1 * zlist[z] * (fx**2))
        H = np.fft.fftshift(H)
        U1=np.fft.fft(np.fft.fftshift(u3x))
        U2 = H * U1
        xdefocus = np.fft.ifftshift(np.fft.ifft(U2))
        xdefocus = np.abs(xdefocus)**2
        xdf = np.vstack((xdf,xdefocus))
        z = z + 1
    plt.figure()
    plt.plot(xdf)
    
    return x1c,x2c,x3c

#def defocus()
    
def display(x1c,x2c,x3c,xdf):
    """displays plots at the focal plane at different wavelengths
    """
    #I NEED to figure out how to initialize truly empty array to avoid this steps
    
    x1c = x1c[1:]
    x2c = x2c[1:]
    x3c = x3c[1:]
    x = 1
    plt.figure()
    while x <= len(x1c):
        plt.plot(x1c[x - 1])
        x += 1
    plt.savefig('plot1.png',format='png')
    
    x = 1
    plt.figure()
    while x <= len(x2c):
        plt.plot(x2c[x - 1])
        x += 1
    plt.savefig('plot2.png',format='png')
    
    x = 1
    plt.figure()
    while x <= len(x3c):
        plt.plot(x3c[x - 1])
        x += 1
    plt.savefig('plot3.png',format='png')
    
    plt.figure()
    y = 1
    p = 0
    zzz = zs - 1
    xdf=np.delete(xdf,0,0)
    dfsum = np.empty([zs-1,M-1])
    while zzz <= len(xdf)+1:
        dfsum = dfsum + xdf[p:zzz]
        p = zzz
        y = y + 1
        zzz = y * (zs - 1)
    plt.imshow(dfsum,aspect='auto')
    
#initialize()
outSamp = np.empty([M-1])
outFoc = np.empty([M-1])
outDef = np.empty([M-1])
lamT=9
out1=np.empty([M-1,lamT])
out2=np.empty([M-1,lamT])
out3=np.empty([M-1,lamT])
'''
c=2.99792458*10**8
lamC=800*10**-9
omega0=2*math.pi*c/lamC
omegar=0.05*10**15
nwavelengths = 2*10**8
frequencylist = np.linspace(omega0-omegar,omega0+omegar,nwavelengths)
i = 0
while i < len(wavelengthlist):
    out1[:,i],out2[:,i],out3[:,i] = prop(wavelengthlist[i])
    i += 1
    
'''         
    
"""
def prop_wave(wavelength):
    prop(wavelength)

if __name__=='__main__':
    pool = Pool()
    pool.map(prop_wave,wavelengthlist)
"""
display(x1c,x2c,x3c,xdf)


