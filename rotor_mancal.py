# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 16:30:12 2020

@author: Lucas
"""
import numpy as np

def rotor(t,h):
  
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
    jg = m/2*((dd/2)**2)
    Tr = 0.1
    wn = np.sqrt(ke/m)
    # dm = 30E-3
    # rm = dm/2
    # lm = 18E-3
    # cr = 100E-6
    # neta = 0.08 
    # F0 = P/2
    
    omega = h[5]
    
    kyy = (-2.44e-17*omega**10 + 4.123e-14*omega**9 - 3.007e-11*omega**8 +
            1.24e-08*omega**7 - 3.18e-06*omega**6 + 0.0005261*omega**5 - 
            0.05641*omega**4 + 3.845*omega**3 - 159.6*omega**2 + 3757*omega + 
            1.105e+05)

    kyz = (7.915e-18*omega**10 - 1.274e-14*omega**9 + 8.724e-12*omega**8 -
            3.294e-09*omega**7 + 7.403e-07*omega**6 - 9.801e-05*omega**5 +
            0.006536*omega**4 + 0.003838*omega**3 - 35.36*omega**2 +
            8223*omega - 2.584e+04)

    kzy = (-1.555e-16*omega**10 + 2.598e-13*omega**9 - 1.87e-10*omega**8 +
            7.589e-08*omega**7 - 1.909e-05*omega**6 + 0.003082*omega**5 -
            0.3198*omega**4 + 20.85*omega**3 - 815*omega**2 +
            1.249e+04*omega - 3.051e+05)

    kzz = (2.47484756e-16*omega**10 - 4.14547720e-13*omega**9 + 2.99180576e-10*omega**8 - 
            1.21728039e-07*omega**7 + 3.06798647e-05*omega**6 - 4.95680070e-03*omega**5 + 
            5.13809382e-01*omega**4 - 3.32787367e+01*omega**3 + 1.27191601e+03*omega**2 - 
            2.59191338e+04*omega + 3.08234066e+05)
    
    cyy = (5.94647821e-17*omega**10 - 9.92150205e-14*omega**9 + 7.12556192e-11*omega**8 - 
            2.88129676e-08*omega**7 + 7.20331783e-06*omega**6 - 1.15105951e-03*omega**5 + 
            1.17454907e-01*omega**4 - 7.42771972*omega**3 + 2.72941946e+02*omega**2 - 
            5.18890960e+03*omega + 5.31875845e+04)    
    
    cyz = (-1.35599625e-16*omega**10 + 2.25971624e-13*omega**9 - 1.62050159e-10*omega**8 +
            6.54031294e-08*omega**7 - 1.63105844e-05*omega**6 + 2.59758374e-03*omega**5 -
            2.63771755e-01*omega**4 + 1.65541148e+01*omega**3 - 6.00162989e+02*omega**2 +
            1.10781530e+04*omega - 8.16874398e+04)
    
    czy = (-1.35599625e-16*omega**10 + 2.25971624e-13*omega**9 - 1.62050159e-10*omega**8 +
            6.54031294e-08*omega**7 - 1.63105844e-05*omega**6 + 2.59758374e-03*omega**5 -
            2.63771755e-01*omega**4 + 1.65541148e+01*omega**3 - 6.00162989e+02*omega**2 +
            1.10781530e+04*omega - 8.16874398e+04)
    
    czz = (6.21354673e-16*omega**10 - 1.03415943e-12*omega**9 + 7.40464964e-10*omega**8 - 
            2.98257250e-07*omega**7 + 7.41873190e-05*omega**6 - 1.17727985e-02*omega**5 + 
            1.18930376*omega**4 - 7.40359150e+01*omega**3 + 2.64563775e+03*omega**2 - 
            4.73091717e+04*omega + 3.24693699e+05)
    
    
    A = 1/jg*(Tr-e*c*h[3]*np.sin(h[2])+e*c*h[4]*np.cos(h[2])-ke*e*h[0]*np.sin(h[2])-ke*e*h[6]*np.sin(h[2])+ke*e*h[1]*np.cos(h[2])-ke*e*h[7]*np.cos(h[2]))
    
    g = np.zeros((8))
    
    g[0] = h[3]
    g[1] = h[4]
    g[2] = h[5]
    g[3] = (1/m*(-c*h[3]-ke*h[0]+ke*h[6]+me*(A*np.sin(h[2])+((h[5])**2)*np.cos(h[2]))))
    g[4] = (1/m*(-P-c*h[4]-ke*h[1]+ke*h[7]+me*(-A*np.cos(h[2])+((h[5])**2)*np.sin(h[2]))))
    g[5] = A
    g[6] = h[6]*(-4*czy*kyz/(4*cyy*czz - 4*cyz*czy) + 2*czz*(ke + 2*kyy)/(4*cyy*czz - 4*cyz*czy)) + h[7]*(-2*czy*(ke + 2*kzz)/(4*cyy*czz - 4*cyz*czy) + 4*czz*kzy/(4*cyy*czz - 4*cyz*czy))-(-2*ke*h[0]*czz/(4*cyy*czz - 4*cyz*czy) + 2*ke*h[1]*cyz/(4*cyy*czz - 4*cyz*czy))
    g[7] = h[6]*(4*cyy*kyz/(4*cyy*czz - 4*cyz*czy) - 2*cyz*(ke + 2*kyy)/(4*cyy*czz - 4*cyz*czy)) + h[7]*(2*cyy*(ke + 2*kzz)/(4*cyy*czz - 4*cyz*czy) - 4*cyz*kzy/(4*cyy*czz - 4*cyz*czy)) - 2*ke*h[0]*czy/(4*cyy*czz - 4*cyz*czy) - 2*ke*h[1]*cyy/(4*cyy*czz - 4*cyz*czy)
    
    # g[6] = (h[6]*((cyz*kzy)/(cyy*czz - cyz*czy) - (czz*(ke + 2*kyy))/(2*(cyy*czz - cyz*czy)))-
    #         h[7]*((czz*kyz)/(cyy*czz - cyz*czy) - (cyz*(ke + 2*kzz))/(2*(cyy*czz - cyz*czy)))+
    #         (ke*h[0]*czz)/(2*(cyy*czz - cyz*czy))-
    #         (ke*h[1]*cyz)/(2*(cyy*czz - cyz*czy))/((cyy*czz)/(cyy*czz - cyz*czy))) 
    
    # g[7] = (-h[6]*((cyy*kzy)/(cyy*czz - cyz*czy) - (czy*(ke + 2*kyy))/(2*(cyy*czz - cyz*czy)))+
    #         h[7]*((czy*kyz)/(cyy*czz - cyz*czy) - (cyy*(5 + 2*kzz))/(2*(cyy*czz - cyz*czy))) -
    #         (ke*h[0]*czy)/(2*(cyy*czz - cyz*czy)) +
    #         (ke*h[1]*cyy)/(2*(cyy*czz - cyz*czy))/((cyy*czz)/(cyy*czz - cyz*czy)))
    

    return np.reshape(g,[8,])








