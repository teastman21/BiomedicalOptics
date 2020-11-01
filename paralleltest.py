# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 13:29:51 2020

parallel test program

@author: teast
"""

import multiprocessing
from multiprocessing import Pool

print(multiprocessing.cpu_count())

def par(x):
    return x*x

if __name__=='__main__':
    with Pool(multiprocessing.cpu_count()-1) as p:
        print(p.map(par,[1,2,3]))
        