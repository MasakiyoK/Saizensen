### Saizensen, Fig. 5.6:
### Drawing the phase diagram of 2-flavor NJL model for nonzero m
### June 2022 / Masakiyo Kitazawa

import numpy as np
from scipy import optimize
import matplotlib as plt

from plt_setting import *

from NJL_funcs import *
m = 4.

def LocalMin( T , mu , Minit , eps=1. ):
    """ Search for local minimum of omega starting from Minit (bracket) """
    return optimize.minimize_scalar( omega , bracket=Minit , args=(m,T,mu,eps) ).x

def GlobalMin( T , mu , Minit1 , Minit2 , eps=1. ):
    """ global minimum; M=0 is always one of minima """
    Mmid1 = ( 2.*Minit1 + Minit2 ) / 3.
    Mmid2 = ( Minit1 + 2.*Minit2 ) / 3.
    min1 = optimize.minimize_scalar( omega , bracket=(Minit1,Mmid1) , args=(m,T,mu,eps) )
    min2 = optimize.minimize_scalar( omega , bracket=(Mmid2,Minit2) , args=(m,T,mu,eps) )
    return min1.x if min1.fun < min2.fun else min2.x

'''
#Check LocalMin/GlobalMin
grid = np.linspace( 0. , 80. , 100 )
for mu in [ 320. , 330. , 340. ] :
    plt.plot( grid , [LocalMin(T,mu,250.) for T in grid] , label=str(mu) , ls="--" )
    plt.plot( grid , [LocalMin(T,mu,150.) for T in grid] , label=str(mu) , ls="-" )
    plt.plot( grid , [GlobalMin(T,mu) for T in grid] , label=str(mu) , ls=":" )
plt.legend()
plt.show()
'''

def Trace1st( T , M1 , M2 , plot=False ):
    """ function for tracing the line of 1st phase transtion from lower T
        M1, M2 : candidates of local minima """
    Mmid  = ( M1 + M2 ) / 2.
    mu = optimize.bisect( lambda mu : GlobalMin(T,mu,M1,M2)-Mmid , 0. , 350. ) # critical mu
    Mmid1 = ( 4.*M1 + M2 ) / 5.
    Mmid2 = ( M1 + 4.*M2 ) / 5.
    M1 = LocalMin( T , mu , (M1,Mmid1) , 1.e-2 )
    M2 = LocalMin( T , mu , (M2,Mmid2) , 1.e-2 )
    print( "T/mu_q/M1/M2 = " , T , " / " , mu , " / " , M1 , " / " , M2 )
    if plot:
        plt.clf()
        Mdif01 = 0.1 * abs( M2 - M1 )
        grid = np.linspace( M1-Mdif01 , M2+Mdif01 , 55 )
        plt.title( "T/mu_q = {:.2f}".format(T) + " / {:.2f}".format(mu) )
        plt.plot( grid , [omega(M,m,T,mu,1.e-1) for M in grid] )
        plt.pause(.5)
    return M1 , M2 , mu
    
#Critical point search
print( "Searching for critical point along 1st tr.:" )
Tlist , mulist , M1list , M2list = [] , [] , [] , []
def AppendTmuM( T , mu , M1 , M2 ):
    Tlist.append(T)
    mulist.append( mu )
    M1list.append( M1 )
    M2list.append( M2 )
    
M1 , M2 , Tdel = 90. , 300. , 4.
for T in np.arange( 0. , 200. , Tdel ):
    M1 , M2 , mu = Trace1st( T , M1 , M2 , plot=True )
    if abs(M2-M1) < 5. : break
    AppendTmuM( T , mu , M1 , M2 )

M1 , M2 , Tdel = M1list[-1] , M2list[-1] , .5
for T in np.arange( Tlist[-1]+Tdel , 200. , Tdel ):
    M1 , M2 , mu = Trace1st( T , M1 , M2 , plot=True )
    if abs(M2-M1) < 1. : break
    AppendTmuM( T , mu , M1 , M2 )

M1 , M2 , Tdel = M1list[-1] , M2list[-1] , .02
for T in np.arange( Tlist[-1]+Tdel , 200. , Tdel ):
    M1 , M2 , mu = Trace1st( T , M1 , M2 , plot=True )
    AppendTmuM( T , mu , M1 , M2 )
    if abs(M2-M1) < .5 : break

plt.clf()
plt.xlabel( "T" )
plt.ylabel( "M" )
plt.plot( Tlist , M1list , marker="D" )
plt.plot( Tlist , M2list , marker="D" )
plt.pause(.5)
plt.show()
plt.clf()

T_crit , mu_crit = Tlist[-1] , mulist[-1]
print( "\ncritical point (T/mu_q) = " , T_crit , " / " , mu_crit )

#Phase diagram
plt.xlabel( "$\mu_q ~\mathrm{[MeV]}$" )
plt.ylabel( "$T ~\mathrm{[MeV]}$" )
plt.xlim( 0. , 363. )
plt.ylim( 0. , 205. )

plt.plot( mulist , Tlist , lw=4 , color="C1" , label="first order" )
plt.plot( [mu_crit] , [T_crit] , marker="o" , ms=10 , color="C1" , label="critical point" , ls="none" )

plt.text( 260. , 185. , r"$m=4~\mathrm{[MeV]}$" )

plt.legend( loc="lower left" )
plt.savefig("phaseDfinite_m.pdf")
