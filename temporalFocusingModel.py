# -*- coding: utf-8 -*-
"""
Temporal Focusing Model

Created on Wed Sep 30 21:51:19 2020


@author: teast
"""


# import necessary packages 
import numpy as np
import math as math
import matplotlib.pyplot as plt


# Define necessary variables
f1 = 0.5        #collimating lens focal length
f2 = 9*10**-3   #objective lens focal length
lamC = 800*10**(-9) #center wavelength
g = 120000      #grating density
thetaDcenter = 0  #diffracted center wavelength 
m = 1           #diffraction order
w = 0.002       #beam width
L = 1.05
M = 2**16
dx = L / M
x = np.arange(-L/2, L / 2 - dx, dx)

#initialize diffraction angle
thetaI = np.arcsin(m * lamC * g - np.sin(thetaDcenter))



def tilt(uin, L, lam, alpha, theta):
    """ported from tilt function, Voelz pg 90
    """
    m = uin.shape[0]
    n = uin.shape[1]
    dx = L/m
    k = 2 * math.pi / lam
    #need to add the small fraction to make end range inclusive
    #solution for a clean fix are write a handmade range function that
    #includes the end range
    x = np.arange((-L/2), (L/2)-dx + .0000000000000001, dx)
    X,Y = np.meshgrid(x, x)
    uout = uin * np.exp(-1j * k * (X * np.cos(theta) + Y*np.sin(theta)) * np.tan(alpha))
    return uout


#initialize beam
u1 = np.exp(-x ** 2 / (2 * w ** 2))

#propagate to the focal plane
L2 = lamC * f1 / dx
dx2 = lamC * f1 / L

#create coordinates at focal plane
xfoc1 = -L2 / 2 + np.arange(0, M - 1) * dx2

#propagate beam with the fourier transform
u2 = 1/(1j * lamC * f1)*(np.fft.ifftshift(np.fft.fft(np.fft.fftshift(u1))))*dx

#propagate to output focal plane
L3 = lamC * f2 / dx2
dx3 = lamC * f2 / L2
 
#create coordinates at output focal plane
xfoc2 = -L3 / 2 + np.arange(0, M-1) * dx3

#propagate beam with fourier transform
u3 = 1/(1j * lamC * f2)*(np.fft.ifftshift(np.fft.fft(np.fft.fftshift(u2))))*dx2 

#u1time = np.fft.ifftshift(np.fft.fft(np.fft.fftshift(u1)))

I=np.abs(u3)**2

plt.figure()
#plt.imshow(I)

print(u3)









