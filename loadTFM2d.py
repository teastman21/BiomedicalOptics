# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 13:51:29 2020

@author: teast
"""

import numpy as np
import matplotlib.pyplot as plt
global nDefocusM
zs = 30

def formLoad():
    global nDefocusM
    nDefocusM = np.load("defocusM.npy")
formLoad()
print(np.shape(nDefocusM))


#dfsum = np.sum(nDefocusM,axis=2,dtype=np.float32)
#plt.imshow(dfsum,aspect='auto')

#xfsum = np.sum(np.abs(np.power(nDefocusM,4)),axis=2,dtype=np.float32)
#xfsum = np.sum(xfsum,axis=1,dtype=np.float32)
lsum = np.sum(nDefocusM,axis=1,dtype=np.float32)
plt.plot(lsum)
#plt.imshow(nDefocusM)
"""
 plt.plot(nDefocusM[15,:,128])
 plt.imshow(nDefocusM[:,:,128],aspect='auto')
"""