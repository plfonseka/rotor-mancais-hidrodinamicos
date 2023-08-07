# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 10:51:12 2020

@author: Lucas
"""
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

gv = 9.8
de = 10E-3
dd = 100E-3
le = 900E-3
ld = 20E-3
me = 1.5E-4
rho = 7850
vol_d = (np.pi*((dd/2)**2))*ld
m = rho*vol_d
e = me/m
E = 200E9
I = np.pi*(de**4)/64
ke = 48*E*I/le**3
c = ke*1E-4
P = m*gv
# jg = m/2*((dd/2)**2)
# Tr = 0.280

dm = 30E-3
rm = dm/2
lm = 18E-3
cr = 100E-6
neta = 0.08 
F0 = P/2

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

for omega in range(1,316):
    
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
    
curve = np.polyfit(W,c22,10)
r_curve = np.reshape(curve,[11])
poly = np.poly1d(r_curve)

print(curve)

new_x = []
new_y = []

for i in range(0,315):
    new_x.append(i+1)
    calc = poly(i+1)
    new_y.append(calc)
    
fig1, ax1 = plt.subplots()  
ax1.plot(W,c11)
ax1.plot(W,c21)
ax1.plot(W,c12)
ax1.plot(W,c22)
# ax1.plot(W,new_y)