### Saizensen, Fig 2.1:
### QCD running coupling (1 loop); Eq.(2.2)
### June 2022 / Masakiyo Kitazawa

import numpy as np
import matplotlib.pyplot as plt
from plt_setting import *

plt.xlabel( r"$\mu/\Lambda_{\rm QCD}$" )
plt.ylabel( r"$\alpha(\mu)=g^2(\mu)/4\pi$" )
plt.xscale("log")
plt.xlim( 0.8 , 800 )
plt.ylim( 0. , 1.3 )
ytics( 0.2 )
plt.xticks( [1   , 3   , 10   , 30   , 100   , 300 ] ,
            ["1" , "3" , "10" , "30" , "100" , "300" ] )

grid_muL = np.exp( np.linspace( 0. , np.log(1000.) , 50 ) )[1:]
alpha = 12.* np.pi / ( 33. - 12. ) / np.log( grid_muL**2 )
plt.plot( grid_muL , alpha )
plt.savefig( "running.pdf" )
