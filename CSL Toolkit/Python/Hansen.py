# updated 20 Nov 2017 by Kerk Phillips

import numpy as np
import matplotlib.pyplot as plt

from LinAppCSL_CSL import LinApp_CSL
from LinAppCSL_SSL import LinApp_SSL
from LinAppCSL_FindSS import LinApp_FindSS
from LinAppCSL_Deriv import LinApp_Deriv
from LinAppCSL_Solve import LinApp_Solve


def Hansen_defs(k,h,z,kp,hp,param):
    A = param[0]
    theta = param[1]
    delta = param[2]

    y = A*(k**theta*(np.exp(z)*hp)**(1-theta))
    i = kp - (1-delta)*k
    c = y - i
    r = theta*y/k
    w = (1-theta)*y/hp

    return y, i, c, r, w


def Hansen_dyn(In, param):
    [kpp, hpp, kp, hp, k, h, zp, z] = In
    hpp =	In[1]
    kp =	In[2]
    hp =	In[3]
    k = 	In[4]
    h = 	In[5]
    zp =	In[6]

    delta = param[2]
    bet = param[3]
    D   = param[4]
    gam = param[5]

    yp, ip, cp, rp, wp = Hansen_defs(kp,hp,zp,kpp,hpp,param)
    y, i, c, r, w = Hansen_defs(k,h,z,kp,hp,param)

    out1 = c**(-gam)*w - D*(1-hp)**(-gam)
    out2 = bet*((c/cp)**gam)*(rp+1-delta)-1
    Out = np.array([out1, out2])
	
    return  Out


# Hansen's model without labor/leisure decision
print("Hansen's model without labor/leisure decision")

#set model parameters
A = 1
theta = .33
delta = .025
bet = .995
gam = 1
rho = .9
sig = .02
D = 2.5
# set up parameter vector to pass to DSGE function file
param = [A, theta, delta, bet, D, gam, rho, sig]

#set numerical parameters
nx = 2
ny = 0
nz = 1
nobs = 250
logX = 1
DO_QZ = False


Zbar = [0]
# find SS numerically
guessXY = [.1, .33]
XYbar = LinApp_FindSS(Hansen_dyn,param,guessXY,Zbar,nx,ny)
print ('XYbar', XYbar)
Xbar = XYbar[0:nx]
Ybar = XYbar[nx:nx+ny]
theta0 = np.append(np.concatenate((Xbar, Xbar, Xbar)),\
            np.concatenate((Zbar, Zbar)) )

NN = rho
#find derivatives and coefficients numerically
AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, WW, TT = \
    LinApp_Deriv(Hansen_dyn,param,theta0,nx,ny,nz,logX)

PP, QQ, UU, RR, SS, VV = \
    LinApp_Solve(AA,BB,CC,DD,FF,GG,HH,JJ,KK,LL,MM,WW,TT,NN,Zbar,0)

print ("PP\n", PP)
print ("QQ\n", QQ)

#find quadratic approximation's steady state\
#[XQtil1, XQtil2] = QuadRoots(.5*HXX,HX-1,H0+sig^2*Hvv/2)

#generate a history of Z's
Z = np.zeros((nobs,nz))
eps = sig*np.random.randn(nobs,nz)
for t in range(1,nobs):
    Z[t,:] = Z[t-1,:].dot(NN) + eps[t,:]

# set starting values and simulate
XYbar = Xbar
X0 = Xbar

XSSL, temp_SSL = LinApp_SSL(X0,Z,XYbar,logX,PP,QQ,UU,[],RR,SS,VV)

XCSL, temp_CSL = LinApp_CSL(Hansen_dyn,param,X0,Z,NN,logX,0,[])

#print ("XSSL", XSSL)
#print ("Difference with CSL",  XSSL - XCSL)

time = range(0, nobs)
plt.figure()
plt.subplot(2,1,1)
plt.plot(time, XSSL[:,0], label='SSL') 
plt.plot(time, XCSL[:,0], label='CSL') 
plt.title('Capital')
plt.legend(loc=2, ncol=1)
plt.xticks([])

plt.subplot(2,1,2)
plt.plot(time, XSSL[:,1], label='SSL') 
plt.plot(time, XCSL[:,1], label='CSL') 
plt.title('Labor')
plt.legend(loc=2, ncol=1)
plt.show()

