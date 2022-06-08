### Saizensen, Fig. 5.3: thermodynamic potential at T>0 , mu=0
### June 2022 / Masakiyo Kitazawa

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
from plt_setting import *

from NJL_funcs import *

# Mininum of omega(M)
M = lambda T , m : abs( optimize.minimize_scalar( omega , bracket=(20.,50.) , args=(m,T,0.) ).x )

### plot
plt.xlabel( r"$T ~\mathrm{[MeV]}$" )
plt.ylabel( r"$M ~\mathrm{[MeV]}$" )
plt.xlim( 0 , 290. )
plt.ylim( 0 , 360. )

m = 4.
grid = np.linspace( 0. , 300. , 50 )
plt.plot( grid , [ M( T,m ) for T in grid ] , ls="--" , label=r"$m=4 ~\mathrm{[MeV]}$" )

m = 0.
Tc = optimize.bisect( lambda T : M(T,m) - 1.e-3 , 150. , 250. )
print( "Tc = " , Tc )

grid = sort( np.concatenate( [ np.linspace( 0. , 300. , 20 ) ,
                               np.linspace( 75. , 168. , 15 ) ,
                               np.linspace( 170. , 192. , 15 ) ,
                               np.linspace( 192.5 , Tc+1.e-3 , 10 ) ] ) )
plt.plot( grid , [ M( T,m ) for T in grid ] , label=r"$m=0 ~\mathrm{[MeV]}$" )

plt.legend( title=r"$\mu_q=0 ~\mathrm{[MeV]}$" )
plt.savefig( "NJL_T_M.pdf" )
