# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 00:58:51 2019

@author: Victor
"""
import random as r
import numpy as np
r.seed(5000)
Nbatch=10000;
Nsample=1250;
x0=np.zeros(Nsample);
y0=np.zeros(Nsample);
xmean=np.zeros(Nbatch);
xave=(2/3)
samplesd=(1/150)
for i in np.arange(0,Nbatch):
    for j in np.arange(1,Nsample):
        k=1;
        while k==1:
            x0[j]=r.random()
            y0[j]=2*r.random()
            if y0[j] < 2*x0[j]:
                k=0
    xmean[i]=np.mean(x0)

sdcount1=0;
sdcount2=0;
sdcount3=0;
for i in np.arange(0,Nbatch):
    
    if float(xmean[i]) < xave+samplesd and float(xmean[i]) > xave-samplesd:
        sdcount1=sdcount1+1

for i in np.arange(0,Nbatch):
    
    if float(xmean[i]) < xave+(2*samplesd) and float(xmean[i]) > xave-(2*samplesd):
        sdcount2=sdcount2+1

for i in np.arange(0,Nbatch):
    
    if float(xmean[i]) < xave+(3*samplesd) and float(xmean[i]) > xave-(3*samplesd):
        sdcount3=sdcount3+1

print(sdcount1/Nbatch,sdcount2/Nbatch,sdcount3/Nbatch)
        