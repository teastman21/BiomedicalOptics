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
        
"""

f1 = 40
d1 = f1
d2 = 20
f2 = -50
d3 = 10
f3 = 30
d4 = 35

import numpy as np
# @TE reminder to use matplotlib to create visual representation of raytracing
import matplotlib.pyplot as plt

#moved this into first function because we can't declare an empty numpy array
#points = np.array([])


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
    pos0 = np.array([int(init_array(angle)[0]),0])
    points = np.array([pos0])
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
    pos1 = np.array([int(d1 + int(position0(angle)[0])), int(array1(angle)[0])])
    points = np.append(points,[pos1],axis=0) 
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
    pos2 = np.array([d2 + int(position1(angle)[0]),int(array2(angle)[0])])
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
    pos3 = np.array([d3 + int(position2(angle)[0]),int(array3(angle)[0])])
    points = np.append(points,[pos3],axis=0)
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
    global points
    pos4 = np.array([d4 + int(position3(angle)[0]), int(array4(angle)[0])])
    points = np.append(points,[pos4],axis=0)
    return pos4

def display():
    """returns a plot of the generated ray tracing
    """
    #@TE i need to make this graphic better... perhaps dynamic?
    x = points[:,0]
    y = points[:,1]
    plt.plot(x,y)
    plt.show()


#cant remember if pi is built into pyth look this up
position4(3.14159/4)
display()




