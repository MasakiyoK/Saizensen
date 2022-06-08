### Saizensen, Fig. 5.5: thermodynamic potential at mu>0 , T=0
### June 2022 / Masakiyo Kitazawa

import matplotlib.pyplot as plt
from plt_setting import *

from NJL_funcs import *

plt.xlabel( r"$M-m ~\mathrm{[MeV]}$" )
plt.ylabel( r"$\omega(M)\times10^{-8} ~\mathrm{[MeV^4]}$" )
plt.xlim( 0. , 380. )
plt.ylim( -1.6 , 0.9 )
xtics( 100. )

mu_list =  [ 335. , 330. , 325. , 310. ]
ls_list =  [ ":"  , "--" , "-." , "-"  ]
col_list = [ "C0" , "C4" , "C3" , "C1" ]

for mu , ls , color in zip( mu_list , ls_list , col_list ) :
    grid = np.linspace( 0. , 380. , 100 )
    plt.plot( grid , [ omega(M,0.,0.,mu)*1.e-8 for M in grid ] , label=r"$\mu_q="+str(round(mu))+" ~\mathrm{[MeV]}$" , ls=ls , color=color )

    minimum = M_T0(mu,0.)
    print( "mu_q / M = " , mu , "/" , minimum )
    plt.plot( [minimum] , [omega(minimum,0.,0.,mu)*1.e-8] , marker="o" , ms=12 , mec=color , color="none" , mew=2 )
    
plt.legend( title=r"$T=0 ~\mathrm{[MeV]}$" )
plt.savefig( "NJL_mu_E.pdf" )
