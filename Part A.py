# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 22:52:36 2019

@author: Victor
"""

import numpy as np
import matplotlib.pyplot as plt
import random as r
r.seed(2014)
tot=1;
scat=0.5;
absp=tot-scat;
q=1;
xb=3.0;
P=(absp/tot);
Nb=10000;
Np=1;
Nc=10;
#x=xb*np.ones(Np);
#sample average of right boundary half range flux and left bound flux
savgrbcurrent=np.zeros(Nb);
savglbflux=np.zeros(Nb);
savgflux1=np.zeros(Nb);
savgflux2=np.zeros(Nb);
savgflux3=np.zeros(Nb);
savgflux4=np.zeros(Nb);
savgflux5=np.zeros(Nb);
savgflux6=np.zeros(Nb);
savgflux7=np.zeros(Nb);
savgflux8=np.zeros(Nb);
savgflux9=np.zeros(Nb);
for k in np.arange(0,Nb):
    x=xb*np.ones(Np);
    rbcount=0;
    lbcount=0;
    count1=0;count2=0;count3=0;count4=0;count5=0;count6=0;count7=0;count8=0;count9=0;
    flux=np.zeros(10);
    for i in np.arange(0,Np):
        in_system=1;
        #sampling initial position (source)
        R1=r.random()
        x[i]=xb*R1
        #sampling initial direction (mu)
        R2=r.random()
        if R2 < 0.5:
            mu=-1
        elif R2 > 0.5:
            mu=1
        #loop will run if particle is still in the system (not absorbed or moved past right boundary)    
        while in_system == 1:
            #sampling distance to collision
            R3=r.random()
            s=-(1/tot)*np.log(1-R3)
            #distance and direction the particle travels before collision
            dx=s*mu
            #cell pathlength counters
            if x[i]<=0.3<x[i]+dx or x[i]>=0.3>x[i]+dx:
                count1=count1+1
            if x[i]<=0.6<x[i]+dx or x[i]>=0.6>x[i]+dx:
                count2=count2+1
            if x[i]<=0.9<x[i]+dx or x[i]>=0.9>x[i]+dx:
                count3=count3+1
            if x[i]<=1.2<x[i]+dx or x[i]>=1.2>x[i]+dx:
                count4=count4+1
            if x[i]<=1.5<x[i]+dx or x[i]>=1.5>x[i]+dx:
                count5=count5+1
            if x[i]<=1.8<x[i]+dx or x[i]>=1.8>x[i]+dx:
                count6=count6+1
            if x[i]<=2.1<x[i]+dx or x[i]>=2.1>x[i]+dx:
                count7=count7+1
            if x[i]<=2.4<x[i]+dx or x[i]>=2.4>x[i]+dx:
                count8=count8+1
            if x[i]<=2.7<x[i]+dx or x[i]>=2.7>x[i]+dx:
                count9=count9+1
            #particle reaches reflective boundary at x=0
            if x[i]+dx <= 0:
                #new position
                x[i]=0
                lbcount=lbcount+1
                mu=-mu
            #particle reaches vaccum boundary at x=xb
            elif x[i]+dx > xb:
                rbcount=rbcount+1
                in_system=0
            #particle remains in system
            else:
                #new position
                x[i]=x[i]+dx
                #sample interaction
                R4=r.random()
                if R4 < P: #absorbed
                    in_system = 0
                elif R4 > P: #scattered
                    #sampling new direction
                    R5=r.random()
                    if R5 < 0.5:
                        mu=-1
                    elif R5 > 0.5:
                        mu=1              
    savgrbcurrent[k]=rbcount*3
    savglbflux[k]=6*lbcount
    savgflux1[k]=3*count1
    savgflux2[k]=3*count2    
    savgflux3[k]=3*count3     
    savgflux4[k]=3*count4    
    savgflux5[k]=3*count5    
    savgflux6[k]=3*count6
    savgflux7[k]=3*count7
    savgflux8[k]=3*count8
    savgflux9[k]=3*count9

flux0=np.mean(savglbflux);
flux1=np.mean(savgflux1);
flux2=np.mean(savgflux2);
flux3=np.mean(savgflux3);
flux4=np.mean(savgflux4);
flux5=np.mean(savgflux5);
flux6=np.mean(savgflux6);
flux7=np.mean(savgflux7);
flux8=np.mean(savgflux8);
flux9=np.mean(savgflux9);
flux10=np.mean(savgrbcurrent);

scalarflux=np.array([flux0,flux1,flux2,flux3,flux4,flux5,flux6,flux7,flux8,flux9,flux10])
xaxis=np.linspace(0,3,11)

plt.plot(xaxis,scalarflux)
plt.title("Scalar flux for Monte-Carlo solution")
plt.xlabel("Position(cm)")
plt.ylabel("Flux Magnitude (n-cm/s)")
plt.show()

sqavgrbcurrent=savgrbcurrent**2;
sqavglbflux=savglbflux**2;
sqavgflux1=savgflux1**2;
sqavgflux2=savgflux2**2;
sqavgflux3=savgflux3**2;
sqavgflux4=savgflux4**2;
sqavgflux5=savgflux5**2;
sqavgflux6=savgflux6**2;
sqavgflux7=savgflux7**2;
sqavgflux8=savgflux8**2;
sqavgflux9=savgflux9**2;

sqflux0=np.mean(sqavglbflux);
sqflux1=np.mean(sqavgflux1);
sqflux2=np.mean(sqavgflux2);
sqflux3=np.mean(sqavgflux3);
sqflux4=np.mean(sqavgflux4);
sqflux5=np.mean(sqavgflux5);
sqflux6=np.mean(sqavgflux6);
sqflux7=np.mean(sqavgflux7);
sqflux8=np.mean(sqavgflux8);
sqflux9=np.mean(sqavgflux9);
sqflux10=np.mean(sqavgrbcurrent);

stdevrbcurrent=np.sqrt((1/(Nb-1)*((sqflux10)-(flux10**2))));
stdevflux0=(1/(Nb-1)*((sqflux0)-(flux0**2)));
stdevflux1=(1/(Nb-1)*((sqflux1)-(flux1**2)));
stdevflux2=(1/(Nb-1)*((sqflux2)-(flux2**2)));
stdevflux3=(1/(Nb-1)*((sqflux3)-(flux3**2)));
stdevflux4=(1/(Nb-1)*((sqflux4)-(flux4**2)));
stdevflux5=(1/(Nb-1)*((sqflux5)-(flux5**2)));
stdevflux6=(1/(Nb-1)*((sqflux6)-(flux6**2)));
stdevflux7=(1/(Nb-1)*((sqflux7)-(flux7**2)));
stdevflux8=(1/(Nb-1)*((sqflux8)-(flux8**2)));
stdevflux9=(1/(Nb-1)*((sqflux9)-(flux9**2)));

relstdevrbcurrent=np.sqrt((1/(Nb-1)*((sqflux10)-(flux10**2))))/0.81451483;
relstdevflux0=(1/(Nb-1)*((sqflux0)-(flux0**2)))/1.71981029;
relstdevflux1=(1/(Nb-1)*((sqflux1)-(flux1**2)))/1.71348235;
relstdevflux2=(1/(Nb-1)*((sqflux2)-(flux2**2)))/1.69421269;
relstdevflux3=(1/(Nb-1)*((sqflux3)-(flux3**2)))/1.66113092;
relstdevflux4=(1/(Nb-1)*((sqflux4)-(flux4**2)))/1.61274277;
relstdevflux5=(1/(Nb-1)*((sqflux5)-(flux5**2)))/1.5468626;
relstdevflux6=(1/(Nb-1)*((sqflux6)-(flux6**2)))/1.46051466;
relstdevflux7=(1/(Nb-1)*((sqflux7)-(flux7**2)))/1.34979871;
relstdevflux8=(1/(Nb-1)*((sqflux8)-(flux8**2)))/1.20971381;
relstdevflux9=(1/(Nb-1)*((sqflux9)-(flux9**2)))/1.03393247;

print(stdevflux0,stdevflux1,stdevflux2,stdevflux3,stdevflux4,stdevflux5,stdevflux6,stdevflux7,stdevflux8,stdevflux9,stdevrbcurrent)

print(relstdevflux0,relstdevflux1,relstdevflux2,relstdevflux3,relstdevflux4,relstdevflux5,relstdevflux6,relstdevflux7,relstdevflux8,relstdevflux9,relstdevrbcurrent)




    