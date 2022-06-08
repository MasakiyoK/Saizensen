### Saizensen, Fig. 8.1:
### Dispersion relations of normal and plasmino modes
### June 2022 / Masakiyo Kitazawa

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from plt_setting import *

def S_Denom( om , p ):
    """ denominator of HTL fermion propagator; See Le Bellac Eq.(6.135) """
    tmp = ( 1. - om/p ) * np.log( (om+p)/(om-p) ) + 2.
    return om - p - 1./2. / p * tmp

E1 = lambda p : optimize.brentq( S_Denom , p+1.e-5 , p+1.1 , args=(p) )
E2 = lambda p : -optimize.brentq( S_Denom , -p-1.e-5 , -p-1. , args=(p) )

### plot
plt.xlabel( r"$p/m_T$" )
plt.ylabel( r"$\omega/m_T$" )
plt.xlim( 0. , 2.1 )
plt.ylim( 0. , 2.55 )
xtics(0.5)

grid = np.linspace( 1.e-5 , 2.1 , 50 )
plt.plot( grid , [ E1(p) for p in grid ] , lw=3 , label="normal" )
plt.plot( grid , [ E2(p) for p in grid ] , lw=3 , ls="--" , label="plasmino" )
plt.plot( grid , grid , ls=":" , lw=2 ) # light cone
plt.fill_between( grid , grid*0 , grid , color='yellow' , alpha=0.2 )
plt.legend()
plt.savefig( "plasmino.pdf" )
