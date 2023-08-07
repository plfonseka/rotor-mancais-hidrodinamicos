# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 18:21:32 2020

@author: Lucas
"""
import numpy as np

def rotor(t,h):
  
    # gv = 9.8
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
    # P = m*gv
    Jg = m/2*((dd/2)**2)
    Tr = 0.004
    
    # dm = 30E-3
    # rm = dm/2
    # lm = 18E-3
    # cr = 100E-6
    # neta = 0.08 
    # F0 = P/2
    
    W = h[5]
    
    kyy =  1e5 *(  0.000000000000009*W**5  - 0.000000000025203*W**4 +  0.000000026498309*W**3  - 0.000012694977000*W**2  + 0.002694675824480*W  + 1.343375683687377);
    kyz =  1e4 *( -0.000000000341610*W**5 +  0.000000296008300*W**4  - 0.000093717642690*W**3 +  0.012970137695722*W**2 -  0.147789396210320*W +  3.843607894856435); 
    kzy =  1e5 *(  0.000000000020692*W**5 -  0.000000019190087*W**4  + 0.000006723597174*W**3 -  0.001101799007177*W**2 +  0.029148376832450*W -  2.511144650787741);
    kzz =  1e5 *(  0.000000000004403*W**5 -  0.000000003041834*W**4  + 0.000000585738878*W**3 +  0.000010958899589*W**2 -  0.012240837879519*W +  1.581158229246111);
   
    cyy =   1e4 *( - 0.000000000000008*W**5 +  0.000000000017949*W**4  - 0.000000022658101*W**3 +  0.000015056407984*W**2 - 0.004637003366861*W  +  1.566556747456347);
    cyz =   1e4 *(   0.000000000000017*W**5  - 0.000000000039915*W**4 +  0.000000050724643*W**3 -  0.000034109679568*W**2 + 0.010781192142140*W  -  1.185639216043245);
    czy =   1e4 *(   0.000000000000017*W**5  - 0.000000000039915*W**4 +  0.000000050724643*W**3 -  0.000034109679568*W**2 + 0.010781192142140*W  -  1.185639216043245);
    czz =   1e4 *(  -0.000000000000055*W**5 +  0.000000000127729*W**4 -  0.000000160671395*W**3 +  0.000106187151324*W**2 - 0.032392726885385*W  +  4.296888513165682);
    
    # Torque aplicado ao conjunto
    # To = 0.300;% Torque aplicado no Eixo
    # 0.280 não passa a zona de ressonância
    # 0.300 passa lentamente a zona de ressonância
    # 0.750 passa rapidamente a zona de ressonância
    
    g = np.zeros((8))
    
    g[0] = h[1]
    g[1] = ((me* h[5]**2 * np.cos(h[4]))/m - h[0]*(ke/m + (e*ke*me*np.sin(h[4])**2)/(Jg*m))-
            h[1]*(c/m + (c*e*me*np.sin(h[4])**2)/(Jg*m)) + (Tr*me*np.sin(h[4]))/(Jg*m)+
            (c* e * me * h[3] * np.cos(h[4])*np.sin(h[4]))/(Jg*m) +
            (e*ke*me*h[2]*np.cos(h[4])*np.sin(h[4]))/(Jg*m)  +  ((me*np.sin(h[4])*(ke*h[7]*e*np.cos(h[4])-
             ke*h[6]*e*np.sin(h[4])))/(Jg*m) - (ke*h[6])/m))
    g[2] = h[3]
    g[3] = ((me* h[5]**2 *np.sin(h[4]))/m - h[2]*(ke/m + (e*ke*me*np.cos(h[4])**2)/(Jg*m))-
            h[3]*(c/m + (c*e*me*np.cos(h[4])**2)/(Jg*m)) - (Tr*me*np.cos(h[4]))/(Jg*m) +
            (c*e*me*h[1]*np.cos(h[4])*np.sin(h[4]))/(Jg*m) + (e*ke*me*h[0]*np.cos(h[4])*np.sin(h[4]))/(Jg*m) -
            ((ke*h[7])/m - (me*np.cos(h[4])*(ke*h[7]*e*np.cos(h[4]) - ke*h[6]*e*np.sin(h[4])))/(Jg*m)))
    g[4] = h[5]
    g[5] = (Tr/Jg + (c*e*h[3]*np.cos(h[4]))/Jg + (e*ke*h[2]*np.cos(h[4]))/Jg - 
            (c*e*h[1]*np.sin(h[4]))/Jg - (e*ke*h[0]*np.sin(h[4]))/Jg + 
            ((ke*h[7]*e*np.cos(h[4]) - ke*h[6]*e*np.sin(h[4]))/Jg))
   
    g[6] = (h[6]*((cyz*kzy)/(cyy*czz - cyz*czy) - (czz*(ke + 2*kyy))/(2*(cyy*czz - cyz*czy)))-
            h[7]*((czz*kyz)/(cyy*czz - cyz*czy) - (cyz*(ke + 2*kzz))/(2*(cyy*czz - cyz*czy)))+
            (ke*h[0]*czz)/(2*(cyy*czz - cyz*czy))-
            (ke*h[3]*cyz)/(2*(cyy*czz - cyz*czy))/((cyy*czz)/(cyy*czz - cyz*czy))) 
    
    g[7] = (-h[6]*((cyy*kzy)/(cyy*czz - cyz*czy) - (czy*(ke + 2*kyy))/(2*(cyy*czz - cyz*czy)))+
            h[7]*((czy*kyz)/(cyy*czz - cyz*czy) - (cyy*(5 + 2*kzz))/(2*(cyy*czz - cyz*czy))) -
            (ke*h[0]*czy)/(2*(cyy*czz - cyz*czy)) +
            (ke*h[3]*cyy)/(2*(cyy*czz - cyz*czy))/((cyy*czz)/(cyy*czz - cyz*czy)))
    
    
    # g[6] = (h[6]*((c12*k21)/(c11*c22 - c12*c21) - (c22*(ke + 2*k11))/(2*(c11*c22 -
    #         c12*c21))) - h[7]*((c22*k12)/(c11*c22 - c12*c21) - 
    #         (c12*(ke + 2*k22))/(2*(c11*c22 - c12*c21))) + (ke*h[0]*c22)/(2*(c11*c22 - c12*c21)) - 
    #         (ke*h[2]*c12)/(2*(c11*c22 - c12*c21)) / ((c11*c22)/(c11*c22 - c12*c21)))
    # g[7] = (-h[6]*((c11*k21)/(c11*c22 - c12*c21) - (c21*(ke + 2*k11))/(2*(c11*c22 - c12*c21))) +
    #         h[7]*((c21*k12)/(c11*c22 - c12*c21) - (c11*(5 + 2*k22))/(2*(c11*c22 - c12*c21))) -
    #         (ke*h[0]*c21)/(2*(c11*c22 - c12*c21)) + (ke*h[2]*c11)/(2*(c11*c22 - c12*c21)) / ((c11*c22)/(c11*c22 - c12*c21)))
    
    return g