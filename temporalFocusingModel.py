# -*- coding: utf-8 -*-
"""
Temporal Focusing Model

Created on Wed Sep 30 21:51:19 2020


@author: Tommy Eastman
"""

"""
Questions for M Durst:
 is diffraction grating tied in properly?
 is this beam right?
 
 
"""

# import necessary packages 
import numpy as np
import math as math
import matplotlib.pyplot as plt


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
    M = 2**10
    
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
    
    #U1X,U1Y = np.meshgrid(u1x,u1y)
    #plot beam after tilt applied
    #I=np.abs(U1X)**2
    plt.figure()
    plt.plot(np.abs(u1x)**2)
    
    plt.figure()
    plt.plot(np.abs(u1y)**2)
    
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
    #U2X,U2Y = np.meshgrid(u2x,u2y)
    #I=np.abs(U2X)**2
    plt.figure()
    plt.plot(np.abs(u2x)**2)
    
    plt.figure()
    plt.plot(np.abs(u2y)**2)
        
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
    #U3X,U3Y = np.meshgrid(u3x,u3y)
    #I=np.abs(U3X)**2
    plt.figure()
    plt.plot(np.abs(u3x)**2)
    
    plt.figure()
    plt.plot(np.abs(u3y)**2)
    
    
    #final coordinates
    #fc = np.meshgrid(u3x,u3y)
        #need to turn fc into x y pairs!
    fc = np.column_stack((u3x,u3y))
    I=np.abs(fc)**2


wavelengthlist = np.linspace(780*10**(-9),820*10**(-9),num=2)
for lamb in wavelengthlist:
    prop(lamb)







