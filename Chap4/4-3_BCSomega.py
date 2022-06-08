### Saizensen, Fig. 4.3: thermodynamic potential for BCS state
### June 2022 / Masakiyo Kitazawa

import numpy as np
from scipy import integrate
from scipy import optimize
import matplotlib.pyplot as plt
from plt_setting import *

N0gv = 0.25
Delta0 = 1./sinh(1./N0gv) # Eq.(4.72)

def omega( Delta , T ):
    """ thermodynamic potential. See Eq.(4.73) """
    def integrand( xi_k ):
        Ek = np.sqrt( xi_k**2 + Delta**2 )
        if( T==0):
            return Ek - xi_k
        ans_D0 = xi_k + 2.* T * log( 1. + exp(-xi_k/T) )
        return Ek + 2.* T * log( 1. + exp(-Ek/T) ) - ans_D0
    integ = integrate.quad( integrand , 0. , 1. )[0]
    return -2.* integ + Delta**2 / N0gv

# minimum of omega at T=0 vs Delta0
print( "minimum of omega: " , 
       Delta0 , 
       optimize.minimize_scalar( omega , bracket=( .5*Delta0 , 2.*Delta0 ) , args=(0.,) ).x )

### plot
plt.xlabel( r"$\Delta/\Delta_0$" )
plt.ylabel( r"$\omega(\Delta)/(N_0 \Delta_0^2)$" )
plt.xlim( 0. , 1.3 )
plt.ylim( -.53 , 0.64 )

def Plot( T , label , ls ):
    grid = np.linspace(0 , 1.5 , 100 )
    result = [ omega( Delta , Delta0*T )/Delta0**2 for Delta in grid*Delta0 ]
    plt.plot( grid , result , label=label , ls=ls )    

vals = ( ( 0.7 , r"$T=0.7\Delta$" , ":" ) ,
         ( 0.6 , r"$T=0.6\Delta$" , (0,(5,1.5,1,1.5,1,1.5)) ) ,
         ( 0.5 , r"$T=0.5\Delta$" , "-." ) ,
         ( 0.3 , r"$T=0.3\Delta$" , "--" ) ,
         ( 0.  , r"$T=0$" , "-" ) )
for T , label , ls in vals:
    Plot( T , label , ls )

plt.legend()
plt.savefig( "BCSomega.pdf" )
