### Saizensen, Fig. 5.2: g dependence of chiral condensate in vacuum 
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
    integ = integrate.quad( integrand , 0. , 1. )[0]
    return -NfNc_pi2 * integ + ( M-m )*( M-m ) / (4.* g )

def MinimumM( m , g ):
    """ Minimum of E """
    return abs( optimize.minimize_scalar( E , bracket=(.05,.3) , args=( m,g*gc ) ).x )

### plot
plt.xlabel( r"$g/g_c$" )
plt.ylabel( r"$(M-m)/\Lambda$" )
plt.xlim( 0.68 , 1.322 )
plt.ylim( 0. , 0.49 )

grid = np.concatenate( [ np.linspace( 0.68 , 0.9  , 20 ) ,
                         np.linspace( 0.9  , 1.0  , 20 ) ,
                         np.linspace( 1.0  , 1.02 , 20 ) ,
                         np.linspace( 1.02 , 1.32 , 20 ) ] )
grid = np.unique( grid )

m_list =  [ 0.004 , 0.002 , 0.  ]
ls_list = [ "-."  , "--"  , "-" ]
for m , ls in zip( m_list , ls_list ) :
    plt.plot( grid , [ MinimumM(m,g)-m for g in grid ] , ls=ls , label=r"$m/\Lambda="+str(m)+"$" )
    
plt.legend()
plt.savefig( "NJL_vac_M.pdf" )
