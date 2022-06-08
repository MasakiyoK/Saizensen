### Saizensen, Fig. 5.3: thermodynamic potential at T>0 , mu=0
### June 2022 / Masakiyo Kitazawa

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
from plt_setting import *

from NJL_funcs import *

plt.xlabel( r"$M-m ~\mathrm{[MeV]}$" )
plt.ylabel( r"$\omega(M)\times10^{-8} ~\mathrm{[MeV^4]}$" )
plt.xlim( 0 , 380. )
xtics( 100. )
plt.ylim( -7.5 , 1.6 )

T_array =  [ 195. , 170. , 120. , 0.   ]
ls_array = [ ":"  , "--" , "-." , "-"  ]
cmap =     [ "C0" , "C4" , "C3" , "C1" ]

for T , ls , color in zip( T_array , ls_array , cmap ) :
    grid = np.linspace( 0. , 380. , 100 )
    plt.plot( grid , [ omega(M,0.,T,0.)*1.e-8 for M in grid ] , label=r"$T="+str(T)+"~\mathrm{[MeV]}$" , ls=ls , color=color )

    minimum = optimize.minimize_scalar( omega , bracket=(50.,200.) , args=(0.,T,0.) )
    print( "T / M = " , T , "/" , minimum.x )
    plt.plot( [minimum.x] , [minimum.fun*1.e-8] , marker="o" , ms=12 , mec=color , color="none" , mew=2 )
    
plt.legend( title=r"$\mu_q=0 ~\mathrm{[MeV]}$" )
plt.savefig( "NJL_T_E.pdf" )
