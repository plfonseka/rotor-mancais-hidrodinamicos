# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 16:13:48 2020

@author: Lucas
"""
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

gv = 9.8
de = 10E-3
le = 800E-3
dd = 90E-3
ld = 15E-3
me = 1E-4
rho = 7850
vol_d = (np.pi*((dd/2)**2))*ld
m = rho*vol_d
e = me/m
E = 200E9
I = np.pi*(de**4)/64
P = m*gv


dm = 30E-3
rm = dm/2
lm = 20E-3
cr = 90E-6
neta = 0.07
F1 = 0.6*P
F2 = 0.4*P
# atribuir F1 p/mancal 1 ou F2 p/mancal 2
F0 = F2


k11 = []
k12 = []
k21 = []
k22 = []

c11 = []
c12 = []
c21 = []
c22 = []

W = []
epsilon = []

for omega in np.arange(0,115):
    
    Fneta = (neta *(lm**3)*omega*rm) / (2*cr**2)
    x0 = 0.7
    
    def S(eps):
        return ((np.pi/2) * (eps/((1-eps**2)**2)) * (np.sqrt(1-eps**2) + (4*eps/np.pi))) - F0/Fneta
    
    feps = fsolve(S,x0)
    
    epsilon.append(feps)

    Ae = 4 / (np.pi**2 + ((16 - np.pi**2)*feps**2))**(3/2)

    aux1 = F0/cr
    
    g11 = (2*np.pi**2 + (16 - np.pi**2) * feps**2) * Ae * aux1
    g12 = (np.pi/4)*((np.pi**2 - 2*np.pi**2 * feps**2 - (16 - np.pi**2)*feps**4 ) / feps*(1 - feps**2)**(1/2)) * Ae * aux1
    g21 = (-np.pi/4)*((np.pi**2 + (32 + np.pi**2)*feps**2 + (32 - 2*np.pi**2)*feps**4)/(feps*(1 - feps**2)**(1/2))) * Ae * aux1
    g22 = ((np.pi**2 + (32 + np.pi**2)*feps**2 + (32 - 2*np.pi**2)*feps**4)/((1 - feps**2)**(1/2))) * Ae * aux1

    k11.append(g11)
    k12.append(g12)
    k21.append(g21)
    k22.append(g22)
        
    aux2 = F0/(cr*omega)
    
    bt11 = (np.pi/2)*(np.sqrt(1 - feps**2) / feps)*(np.pi**2 + (2*np.pi**2 - 16)*feps**2 ) * Ae * aux2
    bt12 = -(2 * np.pi**2 + (4*np.pi**2 - 32)*feps**2 ) * Ae * aux2
    bt21 = bt12
    bt22 = (np.pi/2)*((np.pi**2 + (48 - 2*np.pi**2)*feps**2 + (np.pi**2 * feps**4)) / (np.sqrt(1 - feps**2) * feps)) * Ae * aux2
        
    c11.append(bt11)
    c12.append(bt12)
    c21.append(bt21)
    c22.append(bt22)
    
    W.append(omega)

fig1, ax1 = plt.subplots()    
ax1.plot(W,k11,label='kxx')
ax1.plot(W,k12,label='kxy')
ax1.plot(W,k21,label='kyx')
ax1.plot(W,k22,label='kyy')
ax1.margins(x=0.03)
plt.legend(loc='upper right')
plt.xlabel('Rotação [rad/s]')
plt.ylabel('Rigidez [N/m]')
plt.title('Coeficientes de Rigidez')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.savefig('coef_rig1.png',dpi=600,bbox_inches = 'tight')

fig2, ax1 = plt.subplots()
ax1.plot(W,c11,'b',label='cxx')
ax1.plot(W,c12,'r',label='cxy = cyx')
ax1.plot(W,c22,'c',label='cyy')
ax1.margins(x=0.03)
plt.legend(loc='upper right')
plt.xlabel('Rotação [rad/s]')
plt.ylabel('Amortecimento [Ns/m]')
plt.title('Coeficientes de Amortecimento')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.savefig('coef_amort1.png',dpi=600,bbox_inches = 'tight')
