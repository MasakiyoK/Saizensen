### Saizensen, Fig. 5.2: vacuum energy density
### June 2022 / Masakiyo Kitazawa

import numpy as np
from scipy import integrate
from scipy import optimize
import matplotlib.pyplot as plt
from plt_setting import *

NfNc_pi2 = 6. / np.pi / np.pi
gc = 1./ NfNc_pi2 # Eq.(5.64)

def E( M , m , g ):
    """ Vacuum energy density: Eq. (5.59) """
    integrand = lambda p : p*p * ( np.sqrt( M*M + p*p ) - np.sqrt( m*m + p*p ) )
    integ = integrate.quad( integrand , 0. , 1. )
    return -NfNc_pi2 * integ[0] + ( M-m )*( M-m ) / (4.* g )

### plot
plt.xlabel( r"$(M-m)/\Lambda$" )
plt.ylabel( r"${\cal E}/\Lambda^4\times10^4$" )
plt.xlim( 0 , .31 )
plt.ylim( -5 , 11.5 )
xtics( .1 )

# Make a plot for coupling g and current mass m
def Plot( g , m , ls , color , color_fill , lw , label ):
    grid = np.linspace( 0. , 0.35 , 100 )
    plt.plot( grid-m , [ E(M,m,g*gc)*1.e4 for M in grid ] , label=label , lw=lw , ls=ls , color=color )
    minimum = optimize.minimize_scalar( E , bracket=(.1,.3) , args=(m,g*gc) )
    plt.plot( [minimum.x-m] , [minimum.fun*1.e4] , marker="o" , ms=10 , mew=2 , mec=color , color=color_fill )

plt.plot( [0] , [0] , ls="none" , label = r"$m=0$" ) # dummy plot to just make a legend "m=0"
Plot( 0.9 , 0. , "-." , "C0" , "C0" , 4 , r"$~~~~~~~~~g/g_c=0.9$" )
Plot( 1.0 , 0. , "--" , "C1" , "C1" , 4 , r"$~~~~~~~~~g/g_c=1.0$" )
Plot( 1.1 , 0. , "-"  , "C2" , "C2" , 4 , r"$~~~~~~~~~g/g_c=1.1$" )

plt.plot( [0] , [0] , ls="none" , label = r"$m/\Lambda=0.002$" ) # dummy plot to just make a legend "m/\Lambda=0.002"
Plot( 0.9 , 0.002 , "-." , "C0" , "none" , 2 , r"$~~~~~~~~~g/g_c=0.9$" )
Plot( 1.0 , 0.002 , "--" , "C1" , "none" , 2 , r"$~~~~~~~~~g/g_c=1.0$" )
Plot( 1.1 , 0.002 , "-"  , "C2" , "none" , 2 , r"$~~~~~~~~~g/g_c=1.1$" )
    
plt.legend( ncol=2 , handletextpad=-1.9 )
plt.savefig( "NJL_vac_E.pdf" )
