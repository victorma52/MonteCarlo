# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 22:35:50 2019

@author: Victor
"""
import numpy as np
import matplotlib.pyplot as plt
tot=1
scat=0.5
absp=tot-scat
q=1
xb=3.0
alpha=np.sqrt(absp*tot)
x=np.linspace(0,3,11)
a=-(tot*q)/(absp*((tot*np.cosh(alpha*xb))+(alpha*np.sinh(alpha*xb))))
phi=a*np.cosh(alpha*x)+(q/absp)
plt.plot(x,phi)
plt.title("Scalar flux for Analytical solution")
plt.xlabel("Position(cm)")
plt.ylabel("Flux Magnitude (n-cm/s)")
plt.show()