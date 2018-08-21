# updated 20 Nov 2017 by Kerk Phillips

import numpy as np
import matplotlib.pyplot as plt

from LinAppCSL_CSL import LinApp_CSL
from LinAppCSL_SSL import LinApp_SSL
from LinAppCSL_FindSS import LinApp_FindSS
from LinAppCSL_Deriv import LinApp_Deriv
from LinAppCSL_Solve import LinApp_Solve


def BrockMirman_dyn(In,param):
    kpp =  In[0]
    kp =   In[1]
    k =    In[2]
    zp =   In[3]
    z =    In[4]

    alf = param[0]
    bet = param[1]

    c  = np.exp(z)*k**alf - kp
    cp = np.exp(zp)*kp**alf - kpp
    rp = alf*np.exp(zp)*kp**(alf-1)
    Out = bet*(c/cp)*rp - 1
	
    return Out


##### Borck & Mirman model #####
print("Borck & Mirman model")

#set model parameters
alf = .35
bet = .98
sig = .02
rho = .95
# set up parameter vector to pass to DSGE function file
param = [alf, bet, sig, rho]

#set numerical parameters
nx = 1
ny = 0
nz = 1
nobs = 250
logX = 0


Zbar = [0]
# find SS numerically
XYbar = LinApp_FindSS(BrockMirman_dyn,param,.1,Zbar,nx,ny)
print ('XYbar', XYbar)
Xbar = XYbar[0:nx]
Ybar = XYbar[nx:nx+ny]
theta0 = np.concatenate((Xbar, Xbar, Xbar, Zbar, Zbar))

NN = rho
#find derivatives and coefficients numerically

AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, WW, TT = \
    LinApp_Deriv(BrockMirman_dyn,param,theta0,nx,ny,nz,logX)

PP, QQ, UU, RR, SS, VV = \
    LinApp_Solve(AA,BB,CC,DD,FF,GG,HH,JJ,KK,LL,MM,WW,TT,NN,Zbar,0)

print ("PP\n", PP)
print ("QQ\n", QQ)

#generate a history of Z's
Z = np.zeros((nobs,nz))
eps = sig*np.random.randn(nobs,nz)
for t in range(1,nobs):
    Z[t,:] = Z[t-1,:].dot(NN) + eps[t,:]

# set starting values and simulate
XYbar = Xbar
X0 = Xbar

#  steady state linarization
empty_vec= np.zeros(0)
empty_mat= np.zeros((0,0))

XSSL, temp_SSL = LinApp_SSL(X0,Z,XYbar,logX,PP,QQ,UU,\
                            [],empty_mat,empty_mat,empty_vec)

#  current state linarization
XCSL, temp_CSL = LinApp_CSL(BrockMirman_dyn,param,X0,Z,NN,logX,0,[])

#  exact solution
Xexact = np.zeros((nobs,nx))
Xexact[0,:] = X0
for t in range(1,nobs):
    Xexact[t,:] = alf*bet*np.exp(Z[t,:])*Xexact[t-1,:]**alf

#print ("Exact", Xexact)
#print ("Difference with CSL", Xexact - XCSL)
#print ("Difference with SSL", Xexact - XSSL)

time = range(0, nobs)
plt.figure()
plt.plot(time, Xexact, label = 'Exact') 
plt.plot(time, XSSL, label = 'SSL') 
plt.plot(time, XCSL, label = 'CSL') 
plt.title('Labor')
plt.legend(loc=2, ncol=1)
plt.show()
plt.show()
