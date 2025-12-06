#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 10:45:29 2022

@author: rvwgarcia
"""

import cbm
import time
import numpy as np
import matplotlib.pyplot as plt

NT = int(1e6)
Kin = 1
Taur = 2
Kappa = 3.65
Ps = 0

X0 = np.array([[0.3100, 0.1211]])

tic = time.time()
X, logdX, times, mle, cumulant = cbm.mf(Kin, Taur, Kappa, Ps, NT, X0, 'Lyapunov')
elapsed = time.time() - tic
print(f'Elapsed time: {elapsed} s')

plt.figure(figsize=(12, 8))
plt.subplot(2,1,1)
plt.plot(X[:,0])
plt.grid()
plt.subplot(2,1,2)
plt.plot(times,logdX)
plt.grid()

if len(X)==NT:
    print('BCs are admissible.')
else:
    print('BCs are inadmissible.')


import pandas as pd
import csv
import numpy as np

directory = '/home/rvwgarcia/SublimeScripts/data (1)/'
filename = 'output_0.dat'

with open(directory+filename, newline = '') as file_obj:
    reader_obj = csv.reader(file_obj, delimiter=' ', quotechar=':')
    
    X0 = next(reader_obj)
    
    for row in reader_obj:
        kappa = row[0]
        kappa = kappa[0:-1]
        NT = row[1]
X0 = X0[1:3]




     for row in reader:
         row[0]  # ID
         row[1]  # title
         row[2]  # body
         row[3]  # tags
         
         
         
data = csv.reader(directory+filename, delimeter = ':')

data = pd.read_csv(directory+filename, sep = ':|\s', skiprows=0, nrows=10)

