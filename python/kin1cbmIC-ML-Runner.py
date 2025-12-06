#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 14:41:44 2022

@author: rvwgarcia
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 11:33:41 2022

@author: rvwgarcia
"""

import random
import numpy as np
import torch
import torch.nn as nn
import matplotlib.pyplot as plt

input_size = 4
hidden_size = 6
model_name = 'model_20220929_190927'

num_samples = 100000
X0 = 0.3100
Y0 = 0.1211
KappaMax = 3.7
KappaMin = 3

run_data = [[X0, Y0, (KappaMax-KappaMin)*random.uniform(0, 1)+KappaMin, random.uniform(0, 1)] for i in range(num_samples)]
run_data = torch.tensor(run_data).float()

class ICClassifier(nn.Module):
    def __init__(self, input_size, hidden_size):
        super(ICClassifier, self).__init__()
        self.h1 = nn.Linear(input_size, hidden_size)
        self.out = nn.Linear(hidden_size,1)
    def forward(self, x):
        x = torch.tanh(self.h1(x))
        x = torch.sigmoid(self.out(x))  #use sigmoid with BCELoss
        return x

model = ICClassifier(input_size, hidden_size)

model.load_state_dict(torch.load(model_name))   #load the model state

y_pred = []
for i in range(num_samples):
    # the model on the data
    output = model(run_data[i,:])
               
    #PREDICTIONS
    pred = torch.round(output)
    y_pred.extend(pred.reshape(-1).tolist())

    
run_data = np.array(run_data)
fig = plt.figure()
plt.scatter(run_data[:,3], run_data[:,2], c=y_pred)
plt.colorbar()
plt.show()

            


