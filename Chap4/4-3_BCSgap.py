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

def Minimum( T ):
    """ minimum of omega """
    ans = optimize.minimize_scalar( omega , bracket=( .5*Delta0 , Delta0 ) , args=(T,) ).x
    #Since omega is an even function and ans tends to be negative, take abs
    return abs(ans) if abs(ans)>1.e-3*Delta0 else -1.e-5

#Since Minimum is defined to be negative(-1.e-5) for T>Tc,
#Tc is obtained as follows:
print( "T_c/Delta_0 = " , optimize.bisect( Minimum , 0. , Delta0 )/Delta0 )

### plot
plt.xlabel( r"$T/\Delta_0$" )
plt.ylabel( r"$\Delta/\Delta_0$" )
plt.xlim( 0. , 0.72 )
plt.ylim( 0. , 1.1 )

grid = sort( np.concatenate( [ np.linspace(0 , .72 , 50 ) , np.linspace( 0.55 , 0.57 , 50 ) ] ) )
plt.plot( grid , [ Minimum(T)/Delta0 for T in grid*Delta0 ] )

plt.savefig( "BCSgap.pdf" )
