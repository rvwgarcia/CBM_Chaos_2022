#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 11:33:41 2022

@author: rvwgarcia
"""

import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import TensorDataset
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt

NT = 1e6
input_size = 4
hidden_size = 6

num_samples = 10000
batch_size = 1000
num_batches = int(num_samples/batch_size)

fixedKappas = 1
kappas = np.array([3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.676])
num_kappas = len(kappas)
kappaStr = ['{:.4f}'.format(kappas[i]) for i in range(num_kappas)]

fig_name = 'models_ensemble_BCs.png'
directory = '/home/rvwgarcia/SublimeScripts/'


models = [0]#[1, 3, 5, 7, 9, 11]
#model_names = ['model_kappa'+kappaStr[m] for m in models]
model_names = ['model_20220929_190927'] # 1 hidden layer
#model_names = ['model_20221004_134638'] # 2 hidden layers
Ps = 0

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
        #x = torch.round(x)
        return x

model = ICClassifier(input_size, hidden_size)

#to load model state:
#saved_model = ICClassifier(input_size, hidden_size)


#TESTING THE MODEL
def test(model, test_loader):#, device
    #model in eval mode skips Dropout etc
    model.eval()
    y_true = []
    y_pred = []
    
    # set the requires_grad flag to false as we are in the test mode
    with torch.no_grad():
        for i in test_loader:
            #LOAD THE DATA IN A BATCH
            data,target = i
            
            # moving the tensors to the configured device
            #data, target = data.to(device), target.to(device)
            
            # the model on the data
            output = model(data.float())
                       
            #PREDICTIONS
            pred = np.round(output)
            target = target.float()
            y_true.extend(target.tolist()) 
            y_pred.extend(pred.reshape(-1).tolist())

    #print(f'Accuracy on test set is {accuracy_score(y_true,y_pred)}' , )
    #print("***********************************************************")
    
    return y_true, y_pred#accuracy_score(y_true,y_pred)

for model_name in model_names:
    model.load_state_dict(torch.load(model_name))
    
    if fixedKappas==1:
        accuracy = np.zeros([num_kappas, num_batches])
        colnames = ['X0', 'Y0', 't', 'mle', 'cumulant']
        fig = plt.figure(figsize=(6, 8))
        for i in range(num_kappas):
            filename = 'ScanICs_FixedKappa' + kappaStr[i] + '00.dat'
            #filename = 'ScanBCs_Thu Sep 29 19:32:23 2022.dat'
        
            for j in range(num_batches):
                test_set = pd.read_csv(directory+filename, sep = ' ', names=colnames, header=None, skiprows=j*batch_size, nrows=batch_size)
                labels = test_set.t
                labels = torch.tensor([1 if x==NT else 0 for x in labels]).reshape((batch_size, 1)).float()
                test_set = test_set.drop(columns=['t'])
                
                #if input_size==2:
                test_set = test_set.drop(columns=['mle'])
                test_set = test_set.drop(columns=['cumulant'])
                
                test_set['kappa'] = kappas[i] * np.ones(([batch_size,1]))
                test_set['Ps'] = Ps * np.ones(([batch_size,1]))
                
                test_set = torch.tensor(test_set.values).float()
                test_data = TensorDataset(test_set, labels)   
                test_loader = torch.utils.data.DataLoader(test_data, batch_size, shuffle=True)#batch_size
                Ytrue, YPred = test(model, test_loader)
                accuracy[i][j] = accuracy_score(Ytrue, YPred)
    
        avg_accuracy = np.mean(accuracy, axis=1)
        #print(avg_accuracy)
        std_accuracy = np.std(accuracy, axis=1)
        #plt.imshow(grid, extent=[1,3.6760,0,1], aspect='auto')
        plt.plot(kappas, avg_accuracy, label=model_name)
        plt.fill_between(kappas, (avg_accuracy-std_accuracy), (avg_accuracy+std_accuracy), alpha=0.3)
        plt.xlabel('$\kappa$')
        plt.ylabel('Model accuracy')
        plt.ylim([0,1])
        plt.rcParams.update({'font.size': 12})
        plt.grid()
        plt.legend()
        #plt.savefig(fig_name)
    elif fixedKappas==0:
        filename = 'ScanBCs_Thu Sep 29 19:32:23 2022.dat'

        colnames = ['X0', 'Y0', 'kappa', 'Ps', 'mle', 'cumulant', 't']
        test_set = pd.read_csv(directory+filename, sep = ' ', names=colnames, header=None, skiprows=0, nrows=num_samples)
        labels = test_set.t
        labels = torch.tensor([1 if x==NT else 0 for x in labels]).reshape((num_samples, 1)).float()
        test_set = test_set.drop(columns=['t'])
        test_set = test_set.drop(columns=['mle'])
        test_set = test_set.drop(columns=['cumulant'])
            
        test_set = torch.tensor(test_set.values).float()
        test_data = TensorDataset(test_set, labels)   
        test_loader = torch.utils.data.DataLoader(test_data, num_samples, shuffle=True)#batch_size
        Ytrue, YPred = test(model, test_loader)
        
        accuracy = np.subtract(np.array(Ytrue).reshape((num_samples, 1)),np.array(YPred).reshape((num_samples, 1)))
        
        test_set['accuracy'] = accuracy
        
        dappa = np.array(test_set)
        fig = plt.figure()
        plt.scatter(dappa[:,2], dappa[:,3], c=accuracy)
        plt.colorbar()
        plt.show()

            
#Weights total of 30:
model.h1.weight.data
model.out.weight.data



