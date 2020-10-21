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

    x1c = np.empty([M-1])
    x2c = np.empty([M-1])
    x3c = np.empty([M-1])
    xdf = np.empty([M-1])


def prop(lamm): 
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
    global x1c
    global x2c
    global x3c
    global xdf
    
    #multiplied these two by a hundred due to CPU memory constraints @ output
    dx = L / M
    dy = L / M
    x = np.arange(-L/2, L / 2 - dx + .000000000000000001, dx)
    y = np.arange(-L/2, L / 2 - dy + .000000000000000001, dy)
    k = 2 * math.pi / lamC
    #lam1 = 820 * 10 ** (-9)
    lam1 = lamm
    
    #determine diffraction angle
    thetaI = np.arcsin(m * lamC * g - np.sin(thetaDcenter))
    thetaD = np.arcsin(m * lam1 * g - np.sin(thetaI))
    
    #initialize beam and apply tilt
    u1x = np.exp(-x ** 2 / (2 * w ** 2))* np.exp(-1j * k * (x * np.sin(thetaD)))
    u1y = np.exp(-y ** 2 / (2 * w ** 2))
    

    #plot beam after tilt applied
    plotU1X = np.abs(u1x)**2

    x1c = np.vstack((x1c,plotU1X))
    
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
    x2c = np.vstack((x2c,plotU2X))

        
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

    x3c = np.vstack((x3c,plotU3X))

    
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
    #plt.figure()
    #plt.plot(xdefocus)
   
    
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
    dfsum = np.empty([29,16383])
    while zzz <= len(xdf)+1:
        dfsum = dfsum + xdf[p:zzz]
        p = zzz
        y = y + 1
        zzz = y * (zs - 1)
    plt.imshow(dfsum,aspect='auto')
    
initialize()
wavelengthlist = np.linspace(780*10**(-9),820*10**(-9),num=5)
for lamb in wavelengthlist:
    prop(lamb)
    
display(x1c,x2c,x3c,xdf)



