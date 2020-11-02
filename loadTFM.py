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


dfsum = np.sum(nDefocusM,axis=2,dtype=np.float32)
plt.imshow(dfsum,aspect='auto')