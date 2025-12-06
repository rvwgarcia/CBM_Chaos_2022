#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 13:07:15 2022

@author: rvwgarcia
"""
#import time
import numpy as np

def mf(Kin, Taur, Kappa, Ps, NT, X0, Lyap_Option = 'Basic'):
    X = np.zeros((NT,Taur))
    logdX = np.zeros((NT,1))
    times = np.zeros((NT,1))
    
    X[0,:] = X0
    
    ProdDF = np.identity(Taur)
    m = 0    
    sum2 = 0
    sum4 = 0
    
    for t in range(NT-1):
        SumX = np.sum(X[t,:])
        
        if Lyap_Option=='Lyapunov': #and (NT-t)<=1e4
            dX = np.random.rand(2,1) / 1e6
            norm0 = np.linalg.norm(dX)
    
            DF = cbmJacobian(Kin, Taur, Kappa, Ps, X[t,:], SumX)
            A = DF @ ProdDF
            v = A @ dX
            normt = np.linalg.norm(v)
            
            if normt != 0 and not np.isinf(normt):
                ProdDF = A
                times[m] = t
                logdX[m] = np.log(normt / norm0) / (m+1)
                m += 1

        temp = (1 - SumX) * (Ps + Kappa * (1-Ps) * X[t,0])
        
        if temp >= 0 and temp <= 1 and SumX <= 1:
            X[t+1,0] = temp
            
            for z in range(1,Taur):
                X[t+1,z] = X[t,z-1]
        else:
            print(f'Model constraints violated after {t} iterations.')
            X = X[0:t,:]
            break
    
    if m > 0:
        mle = logdX[m-1]
        
        for n in range(m):
            temp2 = logdX[n] - mle
            sum2 += temp2**2
            sum4 += temp2**4
            
        cumulant = sum4/m - 3 * (sum2/m)**2
        
    if Lyap_Option=='Lyapunov':
        return X, logdX[0:m], times[0:m], mle[0], cumulant[0]
    else:
        return X

def cbmJacobian(Kin, Taur, Kappa, Ps, Xt, SumX):   #for calculating the Jacobian
    DF = np.zeros((Taur,Taur))

    if Kin == 1:
        DF[0,0] = (1 - SumX) * Kappa * (1 - Ps) - (Ps + Kappa * (1-Ps) * Xt[0])
        
        for z in range(1,Taur):
            DF[0,z] = -(Ps + Kappa * (1-Ps) * Xt[0])
            DF[z,z-1] = 1
    elif Kin == 2:
        DF = np.zeros(Taur)

    return DF