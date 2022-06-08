### Saizensen, Fig. 4.2: Fermi distribution function for BCS state
### June 2022 / Masakiyo Kitazawa

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from plt_setting import *

def FermiBCS( e , D ):
    """ Distribution func. of BCS state Eq.(4.52) """
    if D==0. :
        return 0. if e>0. else 1.
    return 1./2. * ( 1. - e/np.sqrt( e*e + D*D ) )

### plot
plt.figure(figsize=(8,4.5))
plt.xlabel( r"$\xi_k$" )
plt.ylabel( r"$\langle S | \hat{c}_{k\sigma}^\dagger \hat{c}_{k\sigma} | S \rangle$" )
plt.xlim( -5 , 5 )
plt.ylim( 0. , 1.05 )

grid = np.linspace( -5. , 5. , 100 )
grid = grid**3 / 25.

plt.plot( grid , [FermiBCS(e,0.) for e in grid] , ls="-" , label=r"$\Delta=0$" )
plt.plot( grid , [FermiBCS(e,1.) for e in grid] , ls="--" , lw=3.5 , label=r"$\Delta=1$")
plt.plot( grid , [FermiBCS(e,2.) for e in grid] , ls=":" , lw=4 , label=r"$\Delta=2$" )

plt.legend()
plt.savefig( "fermi_distBCS.pdf" )
