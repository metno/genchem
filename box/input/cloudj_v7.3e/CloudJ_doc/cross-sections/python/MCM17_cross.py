#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 11:08:19 2023

@author: willemvc
"""

import numpy as np

a = np.loadtxt('/home/willemvc/Desktop/projects/CloudJ/CloudJ_7.3e/cross-sections/MCM17/cross',skiprows=9)

#%%

print(a[:,0])
print(a[:,-1])

#%%

wavelengths = np.arange(202,365.1)
vals = np.interp(wavelengths,a[:,0],a[:,-1])

#%%

b = np.loadtxt('/home/willemvc/Desktop/projects/CloudJ/CloudJ_7.3e/cross-sections/MCM17/qyield',skiprows=7)
print(b[:,0])

vals_q = np.interp(wavelengths,b[:,0],b[:,-1])

#%%

cross = vals * vals_q

#%%
for index, item in enumerate(cross):
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