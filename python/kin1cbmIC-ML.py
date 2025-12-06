#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 13:21:59 2022

@author: rvwgarcia
"""

import pandas as pd
#import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import TensorDataset
#from torch.utils.tensorboard import SummaryWriter
from datetime import datetime
#from sklearn.metrics import accuracy_score

fixedKappa = 0

NT = 1e6
Num_Epochs = 1000
batch_size = 200

input_size = 4
hidden_size = 6

directory = '/home/rvwgarcia/SublimeScripts/'

#1) Separate inputs and targets: First recast data as floats, e.g. ".float()"

if fixedKappa==1:
    kappaStr = '3.1000'
    filename = 'ScanICs_FixedKappa' + kappaStr + '00.dat'
    colnames = ['X0', 'Y0', 't', 'mle', 'cumulant']
elif fixedKappa==0:
    filename = 'ScanBCs_Thu Sep 29 18:34:03 2022.dat'
    colnames = ['X0', 'Y0', 'kappa', 'Ps', 'mle', 'cumulant', 't']
    
training_set = pd.read_csv(directory+filename, sep = ' ', names=colnames, header=None, skiprows=0, nrows=8e3)
validation_set = pd.read_csv(directory+filename, sep = ' ', names=colnames, header=None, skiprows=int(8e3), nrows=2e3)
    
#    
targets = training_set.t
val_targets = validation_set.t
    
targets = torch.tensor([1 if x==NT else 0 for x in targets]).reshape((8000, 1)).float()
val_targets = torch.tensor([1 if x==NT else 0 for x in val_targets]).reshape((2000, 1)).float()
    
training_set = training_set.drop(columns=['t'])
validation_set = validation_set.drop(columns=['t'])

if input_size==2 or fixedKappa==0:
    training_set = training_set.drop(columns=['mle'])
    validation_set = validation_set.drop(columns=['mle'])
    
    training_set = training_set.drop(columns=['cumulant'])
    validation_set = validation_set.drop(columns=['cumulant'])

training_set = torch.tensor(training_set.values).float()
validation_set = torch.tensor(validation_set.values).float()

#2) Create data tensors using TensorDataset:

dataset = TensorDataset(training_set, targets)
val_data = TensorDataset(validation_set, val_targets)

#3) Separate into batches using DataLoader:
    
training_loader = torch.utils.data.DataLoader(dataset, batch_size, shuffle=True)
validation_loader = torch.utils.data.DataLoader(val_data, batch_size, shuffle=True)

#4) Define model:

class ICClassifier(nn.Module):
    def __init__(self, input_size, hidden_size):
        super(ICClassifier, self).__init__()
        self.h1 = nn.Linear(input_size, hidden_size)
        #self.h2 = nn.Linear(hidden_size, hidden_size)
        self.out = nn.Linear(hidden_size,1)
    def forward(self, x):
        x = torch.tanh(self.h1(x))
        #x = torch.tanh(self.h2(x))
        x = torch.sigmoid(self.out(x))  #use sigmoid with BCELoss
        #x = np.round(x)
        return x

model = ICClassifier(input_size, hidden_size)

#loss_fn = torch.nn.MSELoss()
loss_fn = torch.nn.BCELoss()

optimizer = torch.optim.SGD(model.parameters(), lr=0.001, momentum=0.9)

def train_one_epoch(epoch_index):#, tb_writer):
    running_loss = 0
    last_loss = 0
    M = len(training_loader)
    
    for i, data in enumerate(training_loader):
        inputs, labels = data
        
        optimizer.zero_grad()
        outputs = model(inputs)
        loss = loss_fn(outputs, labels)
        loss.backward()     #compute gradients
        optimizer.step()    #update weights
        
        running_loss += loss.item()
        
        if i % M == M-1:
            last_loss = running_loss / M
            #print(f'Batch {i+1} loss: {last_loss}')
    return last_loss

timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
#writer = SummaryWriter('runs/ic_trainer_{}'.format(timestamp))
epoch_number = 0
best_vloss = 1e6

for epoch in range(Num_Epochs):
    # Make sure gradient tracking is on, and do a pass over the data
    model.train(True)
    avg_loss = train_one_epoch(epoch_number)#, writer)

    # We don't need gradients on to do reporting
    model.train(False)

    running_vloss = 0.0
    for i, vdata in enumerate(validation_loader):
        vinputs, vlabels = vdata
        voutputs = model(vinputs)
        vloss = loss_fn(voutputs, vlabels)
        running_vloss += vloss

    avg_vloss = running_vloss / (i + 1)
    
    if (epoch_number + 1) % 100 == 0:
        print(f'EPOCH {epoch_number + 1}: avg train loss = {avg_loss}, avg valid. loss = {avg_vloss}')

    # Log the running loss averaged per batch
    # for both training and validation
    #writer.add_scalars('Training vs. Validation Loss',
    #                { 'Training' : avg_loss, 'Validation' : avg_vloss },
    #                epoch_number + 1)
    #writer.flush()
    
    # Track best performance, and save the model's state
    if avg_vloss < best_vloss:
        best_vloss = avg_vloss
        #model_path = 'model_{}_{}'.format(timestamp, epoch_number)
        #torch.save(model.state_dict(), model_path)
    #model_path = 'model_{}'.format(timestamp)
    
    if fixedKappa==1:
        model_path = 'model_kappa{}_{}'.format(kappaStr,timestamp)
    elif fixedKappa==0:
        model_path = 'model_{}'.format(timestamp)
    torch.save(model.state_dict(), model_path)
    epoch_number += 1

