# -*- coding: utf-8 -*-
"""
Port of Functions from Computational Fourier Optics


Created on Sat Sep 19 21:20:07 2020

@author: Tommy Eastman


"""
import numpy as np
import math as math

x = np.array([[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6],[1,2,3],[4,5,6]])

def propTF(u1,L,lam,z):
    """function ported from 5.1
    takes source field u1 and produces observation field u2
    """
    
    m = u1.shape[0]
    n = u1.shape[1]
    print(m)
    print(n)
    #dx will likely produce a float
    dx = L//m
    print(dx)
    k = 2 * math.pi // lam
    # but the range function cant accept a float?
    print(-1//(2*dx))
    print(1//L)
    print((2*dx)-1//L)
    fx = list(range(-1 // (2 * dx),1//L,(2*dx)-1//L))
    print(fx)