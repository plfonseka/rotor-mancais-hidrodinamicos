# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 20:52:41 2020

@author: Lucas
"""

import numpy as np
import matplotlib.pyplot as plt

cr = 0.0001

r = [2.12591e-05,2.15236e-05,2.10779e-05,
     2.04524e-05,1.98116e-05,1.86368e-05,1.82414e-05,1.71396e-05,
     1.68855e-05,1.60723e-05,1.41368e-05,1.33046e-05,
     1.10025e-05,
     1.05679e-05,8.51674e-06]

r1 = np.array(r)

the = [-1.24711,-1.07487,-0.946036,
       -0.81971,-0.748042,-0.665319,-0.593012,-0.532804,
       -0.530933,-0.518008,-0.475821,-0.405919,
       -0.389015,
       -0.398113,-0.390835]

the1 = np.array(the)

yr = []
zr = []

for tt in np.arange(0,2*np.pi+0.1,0.1):
    y_ = cr*np.sin(tt)
    z_ = cr*np.cos(tt)
    yr.append(y_)
    zr.append(z_)
    
yr1 = np.array(yr)
zr1 = np.array(zr)

def cart2pol(x,y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

rho,theta = cart2pol(yr1,zr1)


# def cart2pol(x,y):
#     rho = np.sqrt(x**2 + y**2)
#     phi = np.arctan2(y, x)
#     return(rho, phi)
  
# # r,the = cart2pol(yr1,zr1)
# rho,theta = cart2pol(r1,the1)

fig = plt.figure()
ax = fig.add_subplot(111, polar=True)
# ax.plot(the,r)
# ax.plot(theta,rho)
ax.plot(the1,r1,ls='none',marker='o',mfc='b',mec='white',ms=15)
ax.set_yticklabels([])
# ax.plot(the1,ys)
# ax.set_theta_zero_location('S')  
plt.savefig('locus_s_folga.png',dpi=600,bbox_inches = 'tight')

