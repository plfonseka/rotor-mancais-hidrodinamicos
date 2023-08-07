# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 11:04:39 2020

@author: Lucas
"""
from sympy import *
from sympy import MatrixSymbol, symbols
from sympy.matrices import Matrix, eye, zeros
from sympy.interactive.printing import init_printing
init_printing(use_unicode=False, wrap_line=False)

To, P, m, th, tht, thp, e, c, cyy, cyz, czy, czz, K, kyy, kyz, kzy, kzz, J, Jg, c, me, ce, Ke, X0P, X1P, X2P, X3P, X4P, X5P, YMp, ZMp, YM, ZM, UY, UZ = symbols('To P m th tht thp e c cyy cyz czy czz K kyy kyz kzy kzz J Jg c me ce Ke X0P X1P X2P X3P X4P X5P YMp ZMp YM ZM UY UZ')

M = Matrix([[m, 0, -me*sin(tht)],[0, m, me*cos(tht)],[0, 0, Jg]]);

Cm = Matrix([[c, 0, -me*cos(tht)*thp],[0, c, -me*sin(tht)*thp],[ce*sin(tht), -ce*cos(tht), 0]])

Km = Matrix([[K, 0, 0],[0, K, 0],[Ke*sin(tht), -Ke*cos(tht), 0]])

F = Matrix([[-K*YM],[-K*ZM + P],[-K*e*sin(tht)*YM + K*e*cos(tht)*ZM + To]])

C1g = Matrix([[2*cyy, 2*cyz],[2*czy, 2*czz]])
K1g = Matrix([[2*kyy+K, 2*kyz],[2*kzy, 2*kzz+K]])
F1g = Matrix([[-K*UY], [-K*UZ]])

MN = M.inv().dot(M)
CN = M.inv().dot(Cm)
KN = M.inv().dot(Km)
FN = M.inv().dot(F)

Der = Matrix([[X0P],[X1P],[X2P],[X3P],[X4P],[X5P]])

C1gD = C1g.inv().dot(C1g)
K1gD = C1g.inv().dot(K1g)
F1gD = C1g.inv().dot(F1g)

V1 = Matrix([[YMp],[ZMp]])
V2 = Matrix([[YM],[ZM]])

aux3 = Matrix(3,3,CN)
aux4 = Matrix(3,3,KN)
C1 = Matrix(2,2,C1gD)
K1 = Matrix(2,2,K1gD)
F1 = Matrix(3,1,FN)
Fg = Matrix(2,1,F1gD)

aux2 = aux4.row_join(aux3)
zero = zeros(3,3)
zeroF = zeros(3,1)

I = eye(3)

aux1 = zero.row_join(I)

Fn = zeroF.col_join(F1)

A = aux1.col_join(aux2)

edo1 = A.dot(Der)
edo2 = C1.dot(V1) 
edo3 = K1.dot(V2)


a = Matrix(3,3,[1,2,1,3,0,6,5,7,9])
b = a.inv()

print(b.dot(a))

'''
Ax
'''
# [X3P,
#  X4P,
#  X5P,
#  K*X0P/m + Ke*X2P*sin(tht)/m + X3P*c/m + X5P*ce*sin(tht)/m,
#  K*X1P/m - Ke*X2P*cos(tht)/m + X4P*c/m - X5P*ce*cos(tht)/m,
#  X2P*(Ke*me*sin(tht)**2/(Jg*m) + Ke*me*cos(tht)**2/(Jg*m)) + X3P*(c*me*sin(tht)/(Jg*m) - me*thp*cos(tht)/Jg) + X4P*(-c*me*cos(tht)/(Jg*m) - me*thp*sin(tht)/Jg) + X5P*(ce*me*sin(tht)**2/(Jg*m) + ce*me*cos(tht)**2/(Jg*m)) + K*X0P*me*sin(tht)/(Jg*m) - K*X1P*me*cos(tht)/(Jg*m)]

'''
F1
'''
# [[0],
#  [0],
#  [0],
#  [-K*YM/m + me*(-K*YM*e*sin(tht) + K*ZM*e*cos(tht) + To)*sin(tht)/(Jg*m)],
#  [(-K*ZM + P)/m - me*(-K*YM*e*sin(tht) + K*ZM*e*cos(tht) + To)*cos(tht)/(Jg*m)],
#  [(-K*YM*e*sin(tht) + K*ZM*e*cos(tht) + To)/Jg]]

'''
C2
'''
# [YMp*(4*cyy*czz/(4*cyy*czz - 4*cyz*czy) - 4*cyz*czy/(4*cyy*czz - 4*cyz*czy)),
#  ZMp*(4*cyy*czz/(4*cyy*czz - 4*cyz*czy) - 4*cyz*czy/(4*cyy*czz - 4*cyz*czy))]

'''
K2
'''
# YM*(-4*czy*kyz/(4*cyy*czz - 4*cyz*czy) + 2*czz*(K + 2*kyy)/(4*cyy*czz - 4*cyz*czy)) + ZM*(-2*czy*(K + 2*kzz)/(4*cyy*czz - 4*cyz*czy) + 4*czz*kzy/(4*cyy*czz - 4*cyz*czy)),
#  YM*(4*cyy*kyz/(4*cyy*czz - 4*cyz*czy) - 2*cyz*(K + 2*kyy)/(4*cyy*czz - 4*cyz*czy)) + ZM*(2*cyy*(K + 2*kzz)/(4*cyy*czz - 4*cyz*czy) - 4*cyz*kzy/(4*cyy*czz - 4*cyz*czy))

'''
F2
'''
# [[-2*K*UY*czz/(4*cyy*czz - 4*cyz*czy) + 2*K*UZ*cyz/(4*cyy*czz - 4*cyz*czy)],
#  [2*K*UY*czy/(4*cyy*czz - 4*cyz*czy) - 2*K*UZ*cyy/(4*cyy*czz - 4*cyz*czy)]]


