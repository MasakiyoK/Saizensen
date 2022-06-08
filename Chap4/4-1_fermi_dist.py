### Saizensen, Fig. 4.1: Fermi distribution function
### June 2022 / Masakiyo Kitazawa

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from plt_setting import *

def Fermi( e , T ):
    """ Fermi distribution func. Eq.(4.27) """
    if T==0. :
        return 0. if e>0. else 1.
    return 1./( exp(e/T) + 1. )

### plot
plt.figure(figsize=(8,4.5))
plt.xlabel( r"$\epsilon-\mu$" )
plt.ylabel( r"$f(\epsilon)$" )
plt.xlim( -5. , 5. )
plt.ylim( 0. , 1.05 )

grid = np.linspace( -5. , 5. , 100 )
grid = grid**3 / 25.

plt.plot( grid , [Fermi(e,0.) for e in grid] , ls="-" , label="$T=0$" )
plt.plot( grid , [Fermi(e,.5) for e in grid] , ls="--" , lw=3.5 , label="$T=0.5$")
plt.plot( grid , [Fermi(e,1.) for e in grid] , ls=":" , lw=4 , label="$T=1$" )

plt.legend()
plt.savefig( "fermi_dist.pdf" , bbox_inches="tight" )
