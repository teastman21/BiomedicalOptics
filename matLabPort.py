# -*- coding: utf-8 -*-
"""
Port of Functions from Computational Fourier Optics


Created on Sat Sep 19 21:20:07 2020

@author: Tommy Eastman


"""
import numpy as np
import math as math
import decimal as d
import cmath as c
import matplotlib.pyplot as plt

x = np.array([[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6]])



def propTF(u1,L,lam,z):
    """function ported from 5.1
    takes source field u1 and produces observation field u2
    """
    m = u1.shape[0]
    n = u1.shape[1]
    dx = L/m
    k = 2 * math.pi / lam
    
    #this combines the initialization of fx and the [FX,FY] step
    fx = np.arange(-1 / (2 * dx),np.abs((2*dx)-1/L),1/L)
    FX,FY=np.meshgrid(fx,fx)
    H = np.exp(1j *math.pi * lam * z * (FX**2 + FY**2))
    H = np.fft.fftshift(H)
    U1=np.fft.fft2(np.fft.fftshift(u1))
    U2 = H * U1
    u2 = np.fft.ifftshift(np.fft.ifft2(U2))
    



def decimal_range(start, stop, step):
    #so return hasn't worked will have to try a while-yield loop
    #need to store start so we can return the final result
    beg = start
    while start <= stop:
        yield float(start)
        start += d.Decimal(step)
    return(list(decimal_range(beg, stop, step)))
#ok this failed yet again

#turns out numpy has a function that did this the whole time :-(
# for i in arange(start,stop,step)
        
 
"""     
lis=[]        
for i in np.arange(-1,100,.66):
    print(i)
    lis.append(i)
"""      

#chapter 6 page 90
#page 208/209
#out=abs(x)<=1/2    
def rect(x):
    out = np.abs(x) <= 1/2
    return out
        
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

plt.figure()
plt.imshow(I1)

u2 = propTF(u1,L1,lam,z)
x2=x1
y2=y1
I2=np.abs(u2)**2

plt.figure()
plt.imshow(I2)

        
    
        
        