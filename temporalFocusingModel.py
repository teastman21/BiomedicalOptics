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
# multiplied these two by a hundred due to CPU memory constraints @ output
dx = 100*L / M
dy = 100* L / M
x = np.arange(-L/2, L / 2 - dx + .000000000000000001, dx)
y = np.arange(-L/2, L / 2 - dy + .000000000000000001, dy)
alpha = 5.0*10**-5 #tilt angle
k = 2 * math.pi / lamC

#determine diffraction angle
thetaI = np.arcsin(m * lamC * g - np.sin(thetaDcenter))


#working with each dimension individually means we dont use this tilt function?
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
    x = np.arange((-L/2), (L/2)-dx + .0000000000000000001, dx)
    X,Y = np.meshgrid(x, x)
    uout = uin * np.exp(-1j * k * (X * np.cos(theta) + Y*np.sin(theta)) * np.tan(alpha))
    return uout


#initialize beam and apply tilt
u1x = np.exp(-x ** 2 / (2 * w ** 2))* np.exp(-1j * k * (x * np.cos(thetaI) + y*np.sin(thetaI)) * np.tan(alpha))
u1y = np.exp(-y ** 2 / (2 * w ** 2))* np.exp(-1j * k * (x * np.cos(thetaI) + y*np.sin(thetaI)) * np.tan(alpha))
U1X,U1Y = np.meshgrid(u1x,u1y)
print(np.shape(u1x))
#plot beam after tilt applied
I=np.abs(U1X)**2
plt.figure()
plt.plot(np.abs(u1x)**2)


if isinstance(u1x, np.ndarray):
    print("yes np")
else:
    print("No np")
#propagate to the focal plane (x)
L2X = lamC * f1 / dx
dx2 = lamC * f1 / L


#y
L2Y = lamC * f1 / dx
dy2 = lamC * f1 / L

#create coordinates at focal plane
xfoc1 = -L2X / 2 + np.arange(0, M - 1) * dx2
yfoc1  = -L2Y / 2 + np.arange(0, M - 1) * dy2

#propagate beam with the fourier transform
u2x = 1/(1j * lamC * f1)*(np.fft.ifftshift(np.fft.fft(np.fft.fftshift(u1x))))*dx
u2y = 1/(1j * lamC * f1)*(np.fft.ifftshift(np.fft.fft(np.fft.fftshift(u1y))))*dy

#plot beam at fourier plane
U2X,U2Y = np.meshgrid(u2x,u2y)
I=np.abs(U2X)**2
plt.figure()
plt.imshow(I)



if isinstance(u2x, np.ndarray):
    print("yes np")
else:
    print("No np")
    
#propagate to output focal plane
L3X = lamC * f2 / dx2
dx3 = lamC * f2 / L2X

L3Y = lamC * f2 / dy2
dy3 = lamC * f2 / L2Y
 
#create coordinates at output focal plane
xfoc2 = -L3X / 2 + np.arange(0, M - 1) * dx3
yfoc2 = -L3Y / 2 + np.arange(0, M - 1) * dy3

#propagate beam with fourier transform
u3x = 1/(1j * lamC * f2)*(np.fft.ifftshift(np.fft.fft(np.fft.fftshift(u2x))))*dx2 
u3y = 1/(1j * lamC * f2)*(np.fft.ifftshift(np.fft.fft(np.fft.fftshift(u2y))))*dy2 

#plot beam at output plane
U3X,U3Y = np.meshgrid(u3x,u3y)
I=np.abs(U3X)**2
plt.figure()
plt.imshow(I)


if isinstance(u3x, np.ndarray):
    print("yes np")
else:
    print("No np")
#final coordinates
#fc = np.meshgrid(u3x,u3y)
    #need to turn fc into x y pairs!
fc = np.column_stack((u3x,u3y))
print(np.shape(fc))
if isinstance(fc, np.ndarray):
    print("yes np")
else:
    print("No np")
    
print(len(fc))
I=np.abs(fc)**2











