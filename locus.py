# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 19:26:05 2020

@author: Lucas
"""

import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter1d

gv = 9.8
de = 10E-3
dd = 100E-3
le = 900E-3
ld = 20E-3
me = 1.5E-4
rho = 7850
vol_d = (np.pi*((dd/2)**2))*ld
m = rho*vol_d
# e = me/m
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

aux1 = []
aux2 = []

for omega in range(1,201,10):

    def orbit(t,h):
                
        Fneta = (neta *(lm**3)*omega*rm) / (2*cr**2)
        x0 = 0.7
        
        def S(eps):
            return ((np.pi/2) * (eps/((1-eps**2)**2)) * (np.sqrt(1-eps**2) + (4*eps/np.pi))) - F0/Fneta
        
        feps = fsolve(S,x0)
        
        Ae = 4 / (np.pi**2 + ((16 - np.pi**2)*feps**2))**(3/2)
    
        aux1 = F0/cr
        
        kyy = (2*np.pi**2 + (16 - np.pi**2) * feps**2) * Ae * aux1
        kyz = (np.pi/4)*((np.pi**2 - 2*np.pi**2 * feps**2 - (16 - np.pi**2)*feps**4 ) / feps*(1 - feps**2)**(1/2)) * Ae * aux1
        kzy = (-np.pi/4)*((np.pi**2 + (32 + np.pi**2)*feps**2 + (32 - 2*np.pi**2)*feps**4)/(feps*(1 - feps**2)**(1/2))) * Ae * aux1
        kzz = ((np.pi**2 + (32 + np.pi**2)*feps**2 + (32 - 2*np.pi**2)*feps**4)/((1 - feps**2)**(1/2))) * Ae * aux1
            
        aux2 = F0/(cr*omega)
        
        cyy = (np.pi/2)*(np.sqrt(1 - feps**2) / feps)*(np.pi**2 + (2*np.pi**2 - 16)*feps**2 ) * Ae * aux2
        cyz = -(2 * np.pi**2 + (4*np.pi**2 - 32)*feps**2 ) * Ae * aux2
        czy = -(2 * np.pi**2 + (4*np.pi**2 - 32)*feps**2 ) * Ae * aux2
        czz = (np.pi/2)*((np.pi**2 + (48 - 2*np.pi**2)*feps**2 + (np.pi**2 * feps**4)) / (np.sqrt(1 - feps**2) * feps)) * Ae * aux2
       
        g = np.zeros((6))
    
        g[0] = h[2]
        g[1] = h[3]
        g[2] = 1/m*(-c*h[2]-ke*h[0]+ke*h[4]+me*((omega**2)*np.cos(omega*t)))
        g[3] = 1/m*(-c*h[3]-ke*h[1]+ke*h[5]+me*((omega**2)*np.sin(omega*t))-P)
        g[4] = h[4]*((cyz*kzy)/(cyy*czz - cyz*czy) - (czz*(ke + 2*kyy))/(2*(cyy*czz - cyz*czy))) - h[5]*((czz*kyz)/(cyy*czz - cyz*czy) - (cyz*(ke + 2*kzz))/(2*(cyy*czz - cyz*czy))) + (ke*h[0]*czz)/(2*(cyy*czz - cyz*czy)) - (ke*h[1]*cyz)/(2*(cyy*czz - cyz*czy))/((cyy*czz)/(cyy*czz - cyz*czy))
        g[5] = -h[4]*((cyy*kzy)/(cyy*czz - cyz*czy) - (czy*(ke + 2*kyy))/(2*(cyy*czz - cyz*czy))) + h[5]*((czy*kyz)/(cyy*czz - cyz*czy) - (cyy*(ke + 2*kzz))/(2*(cyy*czz - cyz*czy))) - (ke*h[0]*czy)/(2*(cyy*czz - cyz*czy)) + (ke*h[1]*cyy)/(2*(cyy*czz - cyz*czy)) /((cyy*czz)/(cyy*czz - cyz*czy))
        
        return g
    
    x0 = np.array([[0],[0],[0]])
    xp0 = np.array([[0],[0],[0]])
    aux = np.block([[x0],[xp0]])
    y0 = np.reshape(aux,(len(aux),))
    
    tin = 0
    tf = 40
    
    sol = solve_ivp(orbit,[tin,tf],y0);
    
    T = sol.t
    Y = sol.y
    
    # w = 150
    
    l= int(len(T)/2)
    # x = Y[0:,l]
    
    a = Y[4,l]
    b = Y[5,l]
    
    # yr = []
    # zr = []
    print(a)
    # for tt in np.arange(0,2*np.pi+0.1,0.1):
    #     y_ = cr*np.sin(tt)
    #     z_ = cr*np.cos(tt)
    #     yr.append(y_)
    #     zr.append(z_)
        
    # yr1 = np.array(yr)
    # zr1 = np.array(zr)
        
    aux1.append(a)
    aux2.append(b)
    
    # ys = gaussian_filter1d(r1, sigma=2)
# print(aux1)
# print(aux2)   
ys = np.array(aux1)
zs = np.array(aux2)
    
def cart2pol(x,y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)
  
# # r,the = cart2pol(yr1,zr1)
r1,the1 = cart2pol(ys,zs)

fig = plt.figure()
ax = fig.add_subplot(111, polar=True)
# ax.plot(the,r)
ax.plot(the1,r1,'bo')
# ax.plot(the1,ys)
# ax.set_theta_zero_location('S')

