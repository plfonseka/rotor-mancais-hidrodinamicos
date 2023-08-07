# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 18:03:53 2020

@author: Lucas
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from rotor_def import rotor

x0 = np.array([[0],[0],[0],[0]])
xp0 = np.array([[0],[0],[0],[0]])
aux = np.block([[x0],[xp0]])
z0 = np.reshape(aux,(len(aux),))

tin = 0
tf = 30

sol = solve_ivp(rotor,[tin,tf],z0,method='RK45');

T = sol.t
Z = sol.y

dof = len(x0)
x = Z[0,:]
xp = Z[5,:]

i = 0

########### PLOT Uy #############
fig1, ax1 = plt.subplots()
ax1.set_xlabel('tempo [s]')
ax1.set_ylabel('Deslocamento [m]')
ax1.plot(T[i:,],Z[6,i:],'b',label='Ym')

ax2 = ax1.twinx()
ax2.set_ylabel('Velocidade de rotação [rad/s]')
ax2.plot(T[i:,],Z[5,i:],'k--',label='\u03A9')

lines_labels = [ax.get_legend_handles_labels() for ax in fig1.axes]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
# plt.legend(loc='upper right')
ax1.legend(lines, labels, loc='lower right')
ax1.margins(x=0.0)
plt.savefig('ym_rapido.png',dpi=600,bbox_inches = 'tight')

########### PLOT Uz ##############
fig2, ax1 = plt.subplots()
ax1.set_xlabel('tempo [s]')
ax1.set_ylabel('Deslocamento [m]')
ax1.plot(T[i:,],Z[7,i:],'b',label='Zm')

ax2 = ax1.twinx()
ax2.set_ylabel('Velocidade de rotação [rad/s]')
ax2.plot(T[i:,],Z[5,i:],'k--',label='\u03A9')

lines_labels = [ax.get_legend_handles_labels() for ax in fig2.axes]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
# plt.legend(loc='upper right')
ax1.legend(lines, labels, loc='lower right')
ax1.margins(x=0.0)
plt.savefig('zm_rapido.png',dpi=600,bbox_inches = 'tight')

