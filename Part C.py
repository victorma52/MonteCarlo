# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 04:13:17 2019

@author: Victor
"""
import random as r
import numpy as np
r.seed(5000)
Nbatch=10000;
Nsample=1250;
x0=np.zeros(Nsample);
y0=np.zeros(Nsample);
wxmean=np.zeros(Nbatch);
wxave=(2/3);
wsamplesd=0.016865;
for i in np.arange(0,Nbatch):
    for j in np.arange(1,Nsample):
        k=1;
        while k==1:
            x0[j]=(r.random())**(1/3)
            y0[j]=2*r.random()
            if y0[j] < 2*x0[j]:
                k=0
    wxmean[i]=np.mean(x0)

sdcount1=0;
sdcount2=0;
sdcount3=0;
for i in np.arange(0,Nbatch):
    
    if float(wxmean[i]) < wxave+wsamplesd and float(wxmean[i]) > wxave-wsamplesd:
        sdcount1=sdcount1+1

for i in np.arange(0,Nbatch):
    
    if float(wxmean[i]) < wxave+(2*wsamplesd) and float(wxmean[i]) > wxave-(2*wsamplesd):
        sdcount2=sdcount2+1

for i in np.arange(0,Nbatch):
    
    if float(wxmean[i]) < wxave+(3*wsamplesd) and float(wxmean[i]) > wxave-(3*wsamplesd):
        sdcount3=sdcount3+1

print(sdcount1/Nbatch,sdcount2/Nbatch,sdcount3/Nbatch)