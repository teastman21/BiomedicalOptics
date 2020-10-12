# -*- coding: utf-8 -*-
"""
Temporal Focusing Model

Created on Wed Sep 30 21:51:19 2020


@author: Tommy Eastman
"""

"""
This branch version is to pass a 2-d np array (x value and wavelength)
so as to remove the need for a for loop
 
 
"""

# import necessary packages 
import numpy as np
import math as math
import matplotlib.pyplot as plt
from multiprocessing import Pool

global x1c
global x2c
global x3c
M = 2**10

def initialize():
    """initialize the np arrays for the coordinate build up
        only run once at beginning!
        """
    global x1c
    global x2c
    global x3c
    x1c = np.array([])
    x2c = np.array([])
    x3c = np.array([])

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
    global x1c
    global x2c
    global x3c
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

    plotU1X = np.abs(u1x)**2
    x1c = np.append(x1c,plotU1X)

    #plt.plot(np.abs(u1y)**2)
    
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
    plotU2X = np.abs(u2x)**2
    x2c = np.append(x2c,plotU2X)

    #plt.plot(np.abs(u2y)**2)
        
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
    plotU3X = np.abs(u3x)**2
    x3c = np.append(x3c,plotU3X)
    #plt.plot(np.abs(u3y)**2)
    
 
#for use with parallel kernels
initialize()
wavelengthlist = np.linspace(780*10**(-9),820*10**(-9),num=3)
for lamb in wavelengthlist:
    prop(lamb)
#lamb = 800*10**(-9) 

def display(x1c,x2c,x3c):
    """displays plots at the focal plane at different wavelengths
    """
    x = 0
    y = 1
    plt.figure()
    while x <= len(x1c):
        plt.plot(x1c[x:y*M])
        x = y * M
        y = y + 1
    
    x = 0
    y = 1
    plt.figure()
    while x <= len(x2c):
        plt.plot(x2c[x:y*M])
        x = y * M
        y = y + 1
        
    x = 0
    y = 1
    plt.figure()
    while x <= len(x3c):
        plt.plot(x3c[x:y*M])
        x = y * M
        y = y + 1

display(x1c,x2c,x3c)

"""
if __name__ == '__main__':
    pool = Pool()
    pool.map(prop, wavelengthlist)
"""  





