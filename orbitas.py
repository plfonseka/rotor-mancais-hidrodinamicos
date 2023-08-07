# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 20:59:58 2020

@author: Lucas
"""
import numpy as np
from scipy.optimize import fsolve

def orbit(t,h):
  
    gv = 9.8
    de = 10E-3
    dd = 100E-3
    le = 900E-3
    ld = 20E-3
    me = 1.5E-4
    rho = 7850
    vol_d = (np.pi*((dd/2)**2))*ld
    m = rho*vol_d
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
    
    # omega = np.sqrt(ke/m)
    # omega = 150
    # omega = 60
    omega = 80
    
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
    
    
    
    
    
    
    
    
    
    
    