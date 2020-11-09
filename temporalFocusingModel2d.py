# -*- coding: utf-8 -*-
"""
Temporal Focusing Model

Created on Wed Sep 30 21:51:19 2020


@author: Tommy Eastman
"""

# import necessary packages 
import numpy as np
import math as math
import matplotlib.pyplot as plt
from multiprocessing import Pool
from datetime import datetime
import multiprocessing
start = datetime.now()



global x1c
global x2c
global x3c
global xdf
global defocusM
global zs

zs = 30 
M = 2**14
lamT=2*4

#initialize necessary arrays
xdf = np.empty([M-1])
outSamp = np.empty([M-1],dtype=np.float32)
outFoc = np.empty([M-1],dtype=np.float32)
outDef = np.empty([M-1],dtype=np.float32)
out1=np.empty([M-1,lamT],dtype=np.float32)
out2=np.empty([M-1,lamT],dtype=np.float32)
out3=np.empty([M-1,lamT],dtype=np.float32)

#frequency information
lambda0 = 800*10**-9
c = 2.99792458*10**8
omega0 = 2 * math.pi * c / lambda0
omegar = 0.05 *10**15
omega = np.linspace(omega0-omegar,omega0+omegar,num=lamT)
lamIn=2*math.pi*c/omega

defocusM = np.empty([zs,lamT],dtype=np.float32)


def prop(lamm,i): 
    # Define necessary variables
    f1 = 0.5        #collimating lens focal length
    f2 = 9*10**-3   #objective lens focal length
    lamC = 800*10**(-9) #center wavelength
    g = 120000      #grating density
    thetaDcenter = 0  #diffracted center wavelength 
    m = 1           #diffraction order
    w = 0.002       #beam width
    L = 1.05
    
    global zs       #defocus slices
    global xdf
    global dx3
    global L3X
    global lam1
    global u3x
    global defocusM
    x1c = np.array([])
    x2c = np.array([])
    x3c = np.array([])
    
    #multiplied these two by a hundred due to CPU memory constraints @ output
    dx = L / M
    #dy = L / M
    x = np.arange(-L/2, L / 2 - dx + .000000000000000001, dx,dtype=np.float32)
    #y = np.arange(-L/2, L / 2 - dy + .000000000000000001, dy)
    k = 2 * math.pi / lamm
    #lam1 = 820 * 10 ** (-9)
    lam1 = lamm
    
    #determine diffraction angle
    thetaI = np.arcsin(m * lamC * g - np.sin(thetaDcenter),dtype=np.float32)
    thetaD = np.arcsin(m * lam1 * g - np.sin(thetaI),dtype=np.float32)
    #initialize beam and apply tilt
    #look at sydneys code to determine variable for 8*10**-9
    #this will tell us what next step is
    #integrate over X, number of points is much smaller
    #parallel processing
    u1x = np.exp(-x ** 2 / (2 * w ** 2))* np.exp(-1j * k * (x * np.sin(thetaD)))*np.exp(-lam1**2/(2*np.abs(8*10**-9)**2))
    #u1y = np.exp(-y ** 2 / (2 * w ** 2))
    

    #plot beam after tilt applied
    plotU1X = np.abs(u1x,dtype=np.float32)**2

    #x1c = np.vstack((x1c,plotU1X))
    x1c = plotU1X
    
    #propagate to the focal plane (x)
    L2X = lam1 * f1 / dx
    dx2 = lam1 * f1 / L
    
    
    #y
    #L2Y = lam1 * f1 / dx
    #dy2 = lam1 * f1 / L
    
    #create coordinates at focal plane
    xfoc1 = -L2X / 2 + np.arange(0, M - 1) * dx2
    #yfoc1  = -L2Y / 2 + np.arange(0, M - 1) * dy2
    
    #propagate beam with the fourier transform
    u2x = 1/(1j * lam1 * f1)*(np.fft.ifftshift(np.fft.fft(np.fft.fftshift(u1x))))*dx
    #u2y = 1/(1j * lam1 * f1)*(np.fft.ifftshift(np.fft.fft(np.fft.fftshift(u1y))))*dy
    #plot beam at fourier plane
  
    plotU2X = np.abs(u2x,dtype=np.float32)**2
    #x2c = np.vstack((x2c,plotU2X))
    x2c = plotU2X
        
    #propagate to output focal plane
    L3X = lam1 * f2 / dx2
    dx3 = lam1 * f2 / L2X
    
    #L3Y = lam1 * f2 / dy2
    #dy3 = lam1 * f2 / L2Y
     
    #create coordinates at output focal plane
    xfoc2 = -L3X / 2 + np.arange(0, M - 1) * dx3
    #yfoc2 = -L3Y / 2 + np.arange(0, M - 1) * dy3
    
    #propagate beam with fourier transform
    u3x = 1/(1j * lam1 * f2)*(np.fft.ifftshift(np.fft.fft(np.fft.fftshift(u2x))))*dx2 
    #u3y = 1/(1j * lam1 * f2)*(np.fft.ifftshift(np.fft.fft(np.fft.fftshift(u2y))))*dy2 
    #plot beam at output plane
    
    plotU3X = np.abs(u3x,dtype=np.float32)**2

    #x3c = np.vstack((x3c,plotU3X))
    x3c = plotU3X
    #defocusCalc()
    defocusM[:,i] = defocusCalc()
    return x1c,x2c,x3c

def defocusCalc():
    #global xdf
    xdf = np.empty([zs],dtype=np.float32)
    fx = (-1 / (2 * dx3) + np.arange(0,M-1,1)*1/L3X)
    #attempt to defocus
    z = 1
    zrange = 1*10**(-3)
    zlist = np.linspace(-zrange,zrange,zs,dtype=np.float32)
    while z < zs:    
        H = np.exp(1j *math.pi * lam1 * zlist[z] * (fx**2))
        H = np.fft.fftshift(H)
        U1=np.fft.fft(np.fft.fftshift(u3x))
        U2 = H * U1
        xdefocus = np.fft.ifftshift(np.fft.ifft(U2))
        xdefocus = np.sum(np.abs(xdefocus)**2)
        xdf[z] = xdefocus
        print(xdf[z])
        z = z + 1
    return xdf
    #plt.figure()
    #plt.plot(xdefocus)

ind = [x for x in range(len(lamIn))]
for i in ind:
    with Pool(multiprocessing.cpu_count()-1) as p:
        out1[:,i],out2[:,i],out3[:,i] = prop(lamIn[i],i)
#xdf=np.delete(xdf,0,0)

#plt.figure()
y = 1
p = 0
zzz = zs - 1
#dfsum = np.sum(defocusM,axis=2,dtype=np.float32)
#plt.imshow(dfsum,aspect='auto')
'''
while zzz < len(defocusM):
    dfsum = dfsum + defocusM[p:zzz]
    p = zzz
    y = y + 1
    zzz = y * (zs - 1)
'''

#plt.imshow(dfsum,aspect='auto')

# Function to save code
def formSave():
    np.save("defocusM.npy",defocusM)
formSave()
print(datetime.now()-start)

"""
def prop_wave(wavelength):
    prop(wavelength)

if __name__=='__main__':
    pool = Pool()
    pool.map(prop_wave,wavelengthlist)

display(x1c,x2c,x3c,xdf)
"""


