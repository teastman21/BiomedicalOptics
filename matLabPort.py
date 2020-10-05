# -*- coding: utf-8 -*-
"""
Port of Functions from Computational Fourier Optics


Created on Sat Sep 19 21:20:07 2020

@author: Tommy Eastman

"""

import numpy as np
import math as math
import decimal as d
import cmath as cmath
import matplotlib.pyplot as plt



def propTF(u1,L,lam,z):
    """function ported from 5.1
    takes source field u1 and produces observation field u2
    """
    m = u1.shape[0]
    n = u1.shape[1]
    dx = L/m
    k = 2 * math.pi / lam
    
    #this combines the initialization of fx and the [FX,FY] step
    fx = np.arange(-1 / (2 * dx), 1/(2*dx)-1/L,1/L)
    FX,FY=np.meshgrid(fx,fx)
    H = np.exp(1j *math.pi * lam * z * (FX**2 + FY**2))
    H = np.fft.fftshift(H)
    U1=np.fft.fft2(np.fft.fftshift(u1))
    U2 = H * U1
    u2 = np.fft.ifftshift(np.fft.ifft2(U2))
    return u2

#chapter 6 page 90
#page 208/209
#out=abs(x)<=1/2    
def rect(x):
    out = np.abs(x) <= 1/2
    return out


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




def focus(uin, L, lam, zf):
    """ported from focus function Voelz pg 94
    """
    m = uin.shape[0]
    n = uin.shape[1]
    dx = L/m
    k = 2 * math.pi / lam
    x = np.arange((-L/2),L/2-dx + .0000001 ,dx)
    X,Y = np.meshgrid(x,x)
    uout = uin * np.exp(1j * k / (2*zf) * (X**2 + Y **2))
    return uout




def propff(u1,L1,lam,z):
    """propFF port from pg 80
    """
    
    m = u1.shape[0]
    n = u1.shape[1]
    dx1 = L1/m
    k = 2*math.pi/lam
    
    L2 = lam*z/dx1
    dx2 = lam*z/L1
    x2 = np.arange(-L2/2, L2/2-dx2, dx2)
    X2,Y2 = np.meshgrid(x2,x2)
    c = 1 / (1j * lam * z) * np.exp(1j * k / (2 * z) * (X2 ** 2 + Y2 **2))
    u2 = c * np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(u1))) * dx1 ** 2
    return u2


#page 66 square beam example
L1 = 0.5
M = 250
dx1 = L1/M
x1 = np.arange((-L1/2),(L1/2 -dx1), dx1)
y1 = x1       
lam = 0.5*10**(-6)
k=2*math.pi/lam
w=0.051
z=2000

X1,Y1 = np.meshgrid(x1,y1)
u1 = rect(X1/(2*w)) * rect(Y1/(2*w))
I1=np.abs(u1**2)
print(np.shape(I1))
print(I1)
plt.figure()
plt.imshow(I1)

"""
deg = math.pi/180
alpha = 5.0 * 10 **-5
theta = 45 *deg
u1 = tilt(u1,L1,lam,alpha,theta)
"""

zf = 2000 
u1 = focus(u1,L1,lam,zf)

u2 = propTF(u1,L1,lam,z)
x2=x1
y2=y1
I2=np.abs(u2)**2

plt.figure()
plt.imshow(I2)

    


        
    
        
        