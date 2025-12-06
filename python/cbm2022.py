"""
cbm2022.py

Runs classic CBM model as described in [Williams-Garcia et al. 2014].

Code uses vectorized operations and modulo arithmetic to minimize runtime.

Updates to come: (1) option to specify desired number of avalanches
                 (2) output in c-webs format
                 (3) implementation of inhibitory nodes (gCBM)
#rvwg2022
"""

#import scipy

#from random import sample
#from random import random
#from random import seed

import random
import numpy as np


N = 20  #number nodes (neurons)
kin = 5 #fixed in-degree
kappa = 1.12   #branching parameter
B = 1.4     #connection weights function coefficient
Taur = 1    #refractory period
Ps = 1e-4   #spontaneous activation probability

NT = 1e5   #fixed number of timesteps (or number of c-webs)

WeightMatrix = np.zeros((N,N))
WeightKernel = np.exp(-B * np.arange(1,kin+1,1))
KernelSum = sum(WeightKernel)

#initialize weight matrix:
for n in range(N):
    Inputs = [*range(N)]
    Inputs.remove(n)    #no self connections
    Inputs = [*map(Inputs.__getitem__,random.sample(range(0,N-1),kin))]
    WeightMatrix[n,Inputs] = kappa * WeightKernel/KernelSum


StateVector = np.zeros(N, dtype='int64')
StateVectorNext = np.zeros(N, dtype='int64')
Activations = np.empty([])

for t in range(NT):
    Actives = [np.where(StateVector == 1)]
    
    if not Actives.size:
        print("empyt")
#Start checking syntax here (8/7/22):
    for n in range(N):
        temp = WeightMatrix[n,Actives != n]

        if len(temp) == 0:
            temp = 0
        else:
            temp = temp-np.random.rand(len(temp))

        StateVector[n,t+1] = ((StateVector[n,t]+1)%(Taur+1))*np.heaviside(np.heaviside(StateVector[n,t],0)+
                                np.heaviside(Ps-random(),0)+sum(np.heaviside(temp,0)),0)

MeanFiringRate = numel(search(StateVector==1))

[i for i,x in enumerate(StateVector) if x==1] 

/(N*NT)

# So yes it can be annoying also... having completions in the comments?

"""
tic
for t=1:NT-1
    Actives = find(StateVector(:,t)==1);
    
    for n=1:N
        temp = WeightMatrix(n,Actives~=n);
        if numel(temp)==0
            temp=0;
        else
            temp = temp-rand(size(temp));
        end
        
        StateVector(n,t+1) = mod(StateVector(n,t)+1,Taur+1)*...
            myStep(myStep(StateVector(n,t))+myStep(Ps-rand)+...
            sum(myStep(temp)));
    end
end
toc

numel(find(StateVector==1))/(N*NT)
"""
