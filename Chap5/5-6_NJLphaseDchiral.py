### Saizensen, Fig. 5.6:
### Drawing the phase diagram of 2-flavor NJL model in the chiral limit
### June 2022 / Masakiyo Kitazawa

import numpy as np
from scipy import optimize
import matplotlib as plt
from plt_setting import *

from NJL_funcs import *

def LocalMin( T , mu , initial , eps=1. ):
    """ Search for local minimum of omega starting from M=initial """
    return abs(optimize.minimize_scalar( omega , bracket=(initial*0.9,initial+1.) , args=(0.,T,mu,eps) ).x)

def GlobalMin( T , mu , initial=250. , eps=1. ):
    """ Global minimum; note that M=0 is always one of extrema at m=0 """
    M = LocalMin( T , mu , initial , eps )
    if M<1.e-3 : return 0.
    return M if omega( M , 0. , T , mu ) < 0. else 0.

'''
#Check LocalMin/GlobalMin
grid = np.linspace( 0. , 80. , 100 )
for mu in [ 300. , 310. , 320. , 330. ] :
    plt.plot( grid , [LocalMin(T,mu,250.) for T in grid] , label=str(mu) , ls="--" )
    plt.plot( grid , [LocalMin(T,mu,0.) for T in grid] , label=str(mu) )
    plt.plot( grid , [GlobalMin(T,mu) for T in grid] , label=str(mu) , ls=":" )
plt.legend()
plt.show()
'''

def Tc( mu , Minitial=250. ):
    """ critical temperature for a given mu """
    if GlobalMin( 0. , mu , Minitial ) == 0. : return 0.
    return optimize.bisect( lambda T : GlobalMin(T,mu,Minitial)-.1 , 0. , 200. )

def mu_c( T , Minitial=250. ):
    """ critical chem. pot. for a given T """
    if GlobalMin( T , 0. , Minitial ) == 0. : return 0.
    return optimize.bisect( lambda mu : GlobalMin(T,mu,Minitial)-.1 , 0. , 350. )

def Trace1st( T , Minitial , plot=False ):
    """ function for tracing phase transtion line from lower T """
    mu = mu_c( T , Minitial )
    M = LocalMin( T , mu , Minitial , 1.e-2 ) #Search minimum starting from Minitial
    print( "T/mu_q/M = " , T , " / " , mu , " / " , M )
    if plot:
        plt.clf()
        grid = np.linspace(0. , M*1.1 , 55 )
        plt.title( "T/mu_q = {:.2f}".format(T) + " / {:.2f}".format(mu) )#str(mu) )
        plt.plot( grid , [omega(M,0.,T,mu,1.e-1) for M in grid] )
        plt.pause(0.5)
    return M , mu
    
#Critical point search
print( "Tracing 1st tr. line:" )
Tlist , mulist , Mlist = [] , [] , []
def AppendTmuM( T , mu , M ):
    Tlist.append(T)
    mulist.append( mu )
    Mlist.append( M )

M , Tdel = 300. , 5.
for T in np.arange( 0. , 200. , Tdel ):
    M , mu = Trace1st( T , .9*M , plot=True )
    if M < 5. : break
    AppendTmuM( T , mu , M )

M , Tdel = Mlist[-1] , .5
for T in np.arange( Tlist[-1]+Tdel , 200. , Tdel ):
    M , mu = Trace1st( T , .9*M , plot=True )
    if M < 2. : break
    AppendTmuM( T , mu , M )
    
M , Tdel = Mlist[-1] , .05
for T in np.arange( Tlist[-1]+Tdel , 200. , Tdel ):
    M , mu = Trace1st( T , .9*M , plot=True )
    AppendTmuM( T , mu , M )
    if M < .5 : break

plt.clf()
plt.xlabel( "T" )
plt.ylabel( "M" )
plt.plot( Tlist , Mlist , marker="D" )
plt.pause(0.5)
plt.show()
plt.clf()

T_crit , mu_crit = Tlist[-1] , mulist[-1]
print( "\ncritical point (T/mu_q) = " , T_crit , " / " , mu_crit )

print( "\nPlotting Phase diagram..." )

grid = np.linspace( 0. , mu_crit , 15 )
plt.plot( grid , [ Tc(mu) for mu in grid ] , ls="--" , label="second order")
plt.plot( mulist , Tlist , lw=4 , label="first order")
plt.plot( [mu_crit] , [T_crit] , marker="o" , ms=10 , ls="none" , color="C1" , label="tri-critical point" )

plt.xlabel( "$\mu_q ~\mathrm{[MeV]}$" )
plt.ylabel( "$T ~\mathrm{[MeV]}$" )
plt.xlim( 0. , 363. )
plt.ylim( 0. , 205. )

plt.text( 140. , 100. , r"$\langle\bar qq\rangle\ne0$" )
plt.text( 275. , 125. , r"$\langle\bar qq\rangle=0$" )
plt.text( 295. , 185. , r"$m=0$" )

plt.legend()
plt.savefig("phaseDchiral.pdf")
