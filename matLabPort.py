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

x = np.array([[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6]])



def propTF(u1,L,lam,z):
    """function ported from 5.1
    takes source field u1 and produces observation field u2
    """
    obs = []
    m = u1.shape[0]
    n = u1.shape[1]
    #dx will likely produce a float
    dx = L/m
    k = 2 * math.pi / lam
    # but the range function cant accept a float?
    print(-1 / (2 * dx))
    print(1 / L)
    print((2*dx)-1/L)
    points = (np.array([[0,0]]))
    
    #this combines the initialization of fx and the [FX,FY] step
    for i in np.arange(-1 / (2 * dx),1/L,abs((2*dx)-1/L)):
        points = np.append(points, [np.array([i,i])],axis=0)
    print(points)
    #need to figure out how to do elementwise operations in np arrays
    h = math.exp(-c.sqrt(-1) * math.pi * lam * z * (points**2).sum(-1))
    print(h)
    #I think it is impossbile to use range with floats in any way, so I will have to hand code a float function
    #Ok this works I am 99% sure just need real numbers to test it with





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
        
        
        

        
        
        
        
        