# -*- coding: utf-8 -*-
"""
Ray Tracing Model

With help from M. Durst's Mathematica Document

By @Tommy Eastman

"""


"""
@TE
Steps forward
    diffraction grating
    this will initalize angle for light
    current framework should work just add diffraction grating functions
    then likely add features like parallel processing
        python multiprocessing package should work well for this
        need to understand difference between synch and asynch
        
        
        allow for multiple wavelengths w for loop
        figure out what the orange line is
        
"""
import numpy as np
# @TE reminder to use matplotlib to create visual representation of raytracing
import matplotlib.pyplot as plt
import math as math
f1 = 0.5
f2 = 0.009
d1 = f1
d2 = f1 + f2
d3 = f2

m = 1
g = 120000
lamC = 800*10**(-9)
thetaDcenter = 0
thetaI = np.arcsin(m * lamC * g - np.sin(thetaDcenter))
print(thetaI)


#moved this into first function because we can't declare an empty numpy array
#points = np.array([])


def diffract(lam):
    """returns diffraction angle
    """
    global m,g,thetaI
    thetaD = np.arcsin(m*lam*g-np.sin(thetaI))
    return thetaD


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
    """@TE sloppy coding I need to clean this up but got stuck in documentation
    but it cant be hard to do!!"""
    
    #clean this up its ugly
    global points
    pos0 = np.array([init_array(angle)[0,0],0])
    print(pos0)
    #points = np.array([pos0])
    #@TE points prints out weird form take a look and fix
    points = [pos0]
    return pos0

def array1(angle):
    """propagate over the distance to d1
    """
    p1 = np.dot(propagation(d1), init_array(angle))
    return p1

def position1(angle):
    """returns x,y coordinates at lens 1 and stores them in the dictionary
    """
    
    """@TE there is a strange space in the output, figure this out"""
    #Learned that numpy adds the strange space for floats in case of a negative sign
    #@TE fix the sloppy code you created while trying to debug here
    global points
    pos1 = np.array([d1 + position0(angle)[0], array1(angle)[0,0]])
    print(points)
    print(pos1)
    points = np.append(points,[pos1],axis=0)
    print(points)
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
    global points
    pos2 = np.array([d2 + position1(angle)[0],array2(angle)[0,0]])
    points = np.append(points,[pos2],axis=0)
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
    global points
    pos3 = np.array([d3 + position2(angle)[0],array3(angle)[0,0]])
    points = np.append(points,[pos3],axis=0)
    print(points)
    return pos3
'''        
def array4(angle):
    """propagates to d4
    """
    
    p4 = np.dot(np.dot(propagation(d4), lens(f3)), array3(angle))
    return p4
    
def position4(angle):
    """
    returns the x,y coordinates at the focal plane and stores them in dic
    """
    global points
    pos4 = np.array([d4 + int(position3(angle)[0]), int(array4(angle)[0])])
    points = np.append(points,[pos4],axis=0)
    return pos4
'''
def display():
    """returns a plot of the generated ray tracing
    """
    #@TE i need to make this graphic better... perhaps dynamic?
    x = points[:,0]
    y = points[:,1]
    plt.plot(x,y)
    plt.show()


wavelengthlist = np.linspace(780*10**(-9),820*10**(-9),num=5)
for lamb in wavelengthlist:
    position3(diffract(lamb))
    display()
    
    
#position3(diffract(820*10**(-9)))
#display()



