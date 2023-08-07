# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 21:11:29 2020

@author: Lucas
"""
import numpy as np
from orbitas import orbit
from scipy.integrate import solve_ivp
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

dm = 30E-3
rm = dm/2
lm = 18E-3
cr = 100E-6
neta = 0.08 
F0 = P/2

x0 = np.array([[0],[0],[0]])
xp0 = np.array([[0],[0],[0]])
aux = np.block([[x0],[xp0]])
y0 = np.reshape(aux,(len(aux),))

tin = 0
tf = 30

sol = solve_ivp(orbit,[tin,tf],y0);

T = sol.t
Y = sol.y

l=1000

x = Y[0:,l:1050]

a = x[4,:]
b = x[5,:]


# yr = []
# zr = []

# for tt in np.arange(0,2*np.pi+0.1,0.1):
#     y_ = cr*np.sin(tt)
#     z_ = cr*np.cos(tt)
#     yr.append(y_)
#     zr.append(z_)
    
# yr1 = np.array(yr)
# zr1 = np.array(zr)
    
def cart2pol(x,y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

# r,the = cart2pol(yr1,zr1)
r1,the1 = cart2pol(a,b)

fig = plt.figure()
ax = fig.add_subplot(111, polar=True)
ax.plot(the1,r1,'k:',label='Órbita do Mancal')
# ax.plot(the,r,'b',label='Folga Radial')
ax.set_yticklabels([])
# fig.legend(loc='best')

plt.savefig('orbit_80w_sf.png',dpi=600,bbox_inches = 'tight')
