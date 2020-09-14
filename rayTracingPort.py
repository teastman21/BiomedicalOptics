# -*- coding: utf-8 -*-
"""
Ray Tracing Model

Ported from M. Durst's Mathematica Document

By @Tommy Eastman

"""


"""
So far I have coordinate values at "critical points" (need to know exactly
what these critical points are) auto-populate into a single dictionary
(is a an array a better idea?) where they can be extracted from later 
"""

f1 = 40
d1 = f1
d2 = 20
f2 = -50
d3 = 10
f3 = 30
d4 = 35

import numpy as np

def propagation(d):
    """applies propagation matrix through a distance d
    """
    
    m = np.array([[1,d],[0,1]])
    return m

def lens(f):
    """applies a matrix for the lens of focal length f
    """
    m = np.array([[1,0],[-1/f,1]])
    return m

def init_array(angle):
    """array containing initial angle and initial position
    """
    m = np.array([[0],[angle]])
    return m

def position0(angle):
    """returns coordinates at initial position
    """
    """@Tommy OK this is sloppy coding need to clean this up but got stuck in documentation
    but it cant be hard to do!!"""
    x = int(init_array(angle)[0])
    c = np.array([x,0])
    return c

def array1(angle):
    """propagate over the distance to d1
    """
    p1 = np.dot(propagation(d1), init_array(angle))
    return p1

def position1(angle):
    """returns x,y coordinates at lens 1 and stores them in the dictionary
    """
    
    """@Tommy there is a strange space in the output, figure this out"""
    
    pos1 = np.array([d1 + int(position0(angle)[0]), int(array1(angle)[0])])
    return pos1

def array2(angle):
    """propagates to d2
    """
    
    p2 = np.dot(np.dot(propagation(d2), lens(f1)), array1(angle))
    return p2
    
def position2(angle):
    """
    returns the x,y coordinates at the second lens and stores them in dic
    """

    pos2 = np.array([d2 + int(position1(angle)[0]),int(array2(angle)[0])])
    return pos2
    
def array3(angle):
    """propagates to d3
    """
    
    p3 = np.dot(np.dot(propagation(d3), lens(f2)), array2(angle))
    return p3
    
def position3(angle):
    """
    returns the x,y coordinates at the third lens and stores them in dic
    """

    pos3 = np.array([d3 + int(position2(angle)[0]),int(array3(angle)[0])])
    return pos3
        
def array4(angle):
    """propagates to d4
    """
    
    p4 = np.dot(np.dot(propagation(d4), lens(f3)), array3(angle))
    return p4
    
def position4(angle):
    """
    returns the x,y coordinates at the focal plane and stores them in dic
    """

    pos4 = np.array([d4 + int(position3(angle)[0]), int(array4(angle)[0])])
    return pos4   


"""running position4 will propagate through the model @Tommy add a way to populate
a dictionary so that we can look back at the critical locations
"""

print(position4(20))




