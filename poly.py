# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 11:26:29 2020

@author: Lucas
"""
import matplotlib.pyplot as plt

k11 = []
k12 = []
k21 = []
k22 = []

c11 = []
c12 = []
c21 = []
c22 = []

w = []

for omega in range(1,200):

    k_11 = (-2.44e-17*omega**10 + 4.123e-14*omega**9 - 3.007e-11*omega**8 +
           1.24e-08*omega**7 - 3.18e-06*omega**6 + 0.0005261*omega**5 - 
           0.05641*omega**4 + 3.845*omega**3 - 159.6*omega**2 + 3757*omega + 
           1.105e+05)

    k_12 = (7.915e-18*omega**10 - 1.274e-14*omega**9 + 8.724e-12*omega**8 -
           3.294e-09*omega**7 + 7.403e-07*omega**6 - 9.801e-05*omega**5 +
           0.006536*omega**4 + 0.003838*omega**3 - 35.36*omega**2 +
           8223*omega - 2.584e+04)

    k_21 = (-1.555e-16*omega**10 + 2.598e-13*omega**9 - 1.87e-10*omega**8 +
           7.589e-08*omega**7 - 1.909e-05*omega**6 + 0.003082*omega**5 -
           0.3198*omega**4 + 20.85*omega**3 - 815*omega**2 +
           1.249e+04*omega - 3.051e+05)

    k_22 = (2.47484756e-16*omega**10 - 4.14547720e-13*omega**9 + 2.99180576e-10*omega**8 - 
            1.21728039e-07*omega**7 + 3.06798647e-05*omega**6 - 4.95680070e-03*omega**5 + 
            5.13809382e-01*omega**4 - 3.32787367e+01*omega**3 + 1.27191601e+03*omega**2 - 
            2.59191338e+04*omega + 3.08234066e+05)
    
    k11.append(k_11)
    k12.append(k_12)
    k21.append(k_21)
    k22.append(k_22)
    
    c_11 = (5.94647821e-17*omega**10 - 9.92150205e-14*omega**9 + 7.12556192e-11*omega**8 - 
            2.88129676e-08*omega**7 + 7.20331783e-06*omega**6 - 1.15105951e-03*omega**5 + 
            1.17454907e-01*omega**4 - 7.42771972*omega**3 + 2.72941946e+02*omega**2 - 
            5.18890960e+03*omega + 5.31875845e+04)    
    
    c_12 = (-1.35599625e-16*omega**10 + 2.25971624e-13*omega**9 - 1.62050159e-10*omega**8 +
            6.54031294e-08*omega**7 - 1.63105844e-05*omega**6 + 2.59758374e-03*omega**5 -
            2.63771755e-01*omega**4 + 1.65541148e+01*omega**3 - 6.00162989e+02*omega**2 +
            1.10781530e+04*omega - 8.16874398e+04)
    
    c_21 = c_12
    
    c_22 = (6.21354673e-16*omega**10 - 1.03415943e-12*omega**9 + 7.40464964e-10*omega**8 - 
            2.98257250e-07*omega**7 + 7.41873190e-05*omega**6 - 1.17727985e-02*omega**5 + 
            1.18930376*omega**4 - 7.40359150e+01*omega**3 + 2.64563775e+03*omega**2 - 
            4.73091717e+04*omega + 3.24693699e+05) 
     
    c11.append(c_11)
    c12.append(c_12)
    c21.append(c_21)
    c22.append(c_22) 
    
    w.append(omega)

fig1, ax1 = plt.subplots()    
ax1.plot(w,k11,label='kxx')
ax1.plot(w,k12,label='kxy')
ax1.plot(w,k21,label='kyx')
ax1.plot(w,k22,label='kyy')
ax1.margins(x=0.03)
plt.legend(loc='upper right')
plt.xlabel('Rotação [rad/s]')
plt.ylabel('Rigidez [N/m]')
plt.title('Rigidez Polinomial')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

plt.savefig('rig_poly.png',dpi=600)

fig2, ax1 = plt.subplots()
ax1.plot(w,c11,'b',label='cxx')
ax1.plot(w,c12,'r',label='cxy = cyx')
ax1.plot(w,c22,'c',label='cyy')
ax1.margins(x=0.03)
plt.legend(loc='upper right')
plt.xlabel('Rotação [rad/s]')
plt.ylabel('Amortecimento [Ns/m]')
plt.title('Amortecimento Polinomial')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

plt.savefig('amort_poly.png',dpi=600)







