### Saizensen, Fig. 5.5: constituent quark mass at mu>0 , T=0
### June 2022 / Masakiyo Kitazawa

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
from plt_setting import *

from NJL_funcs import *

plt.xlabel( r"$\mu_q ~\mathrm{[MeV]}$" )
plt.ylabel( r"$M ~\mathrm{[MeV]}$" )
plt.xlim( 0. , 420. )
plt.ylim( 0. , 345. )

def Plot( m , ls ):
    mu_c = optimize.bisect( lambda mu : M_T0(mu,m)-150. , 300. , 400. )
    print( "Critical mu_q = " , mu_c )
    grid = np.sort( np.concatenate( [ np.linspace(0 , 290. , 10 ) ,
                                      np.linspace( 300. , mu_c-1. , 15 ) ,
                                      np.linspace( mu_c+1. , 420. , 15 ) ,
                                      np.linspace( mu_c-0.01 , mu_c+0.01 , 2 ) ] ) )
    plt.plot( grid , [ M_T0(mu,m) for mu in grid ] , ls=ls , label=r"$m="+str(round(m))+"~\mathrm{[MeV]}$" )
    
Plot( 4. , "--" )
Plot( 0. , "-" )

plt.legend( title=r"$T=0 ~\mathrm{[MeV]}$" )
plt.savefig( "NJL_mu_M.pdf" )
