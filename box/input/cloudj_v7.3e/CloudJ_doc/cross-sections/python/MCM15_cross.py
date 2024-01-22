#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 11:08:19 2023

@author: willemvc
"""

import numpy as np

a = np.loadtxt('/home/willemvc/Desktop/projects/CloudJ/CloudJ_7.3e/cross-sections/MCM15/cross',skiprows=18)

#%%

print(a[:,0])
print(a[:,-1])

#%%

wavelengths = np.arange(202,364.1)
vals = np.interp(wavelengths,a[:,0],a[:,-1])

#%%
for index, item in enumerate(vals):
    print(item)
    # print(int(item))
    # print(int(wavelengths[index]))
    # print(a[index,0], item)

#%%

a = np.loadtxt('/home/willemvc/Desktop/projects/CloudJ/CloudJ_7.3e/cross-sections/GLYOX/cross',skiprows=13)

print(a[:,0])
print(a[:,-1])

for index, item in enumerate(a[:211,-1]):
    print(item)
    # print(int(item))
    
#%%
print(460 - 250)
#%%