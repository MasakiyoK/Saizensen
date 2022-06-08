### Saizensen, Fig. 5.4:
### lattice result of chiral condensate (and Polyakov loop)
### June 2022 / Masakiyo Kitazawa

import numpy as np
import matplotlib.pyplot as plt
from plt_setting import *

# lattice results of chiral condensate and Polyakov loop
# from JHEP 1009, 073 (2010); Table 3
# https://arxiv.org/abs/1005.3508
dat = np.loadtxt( "chiral_pol.dat" )

### plot
plt.xlabel( r"$T \,\, \mathrm{[MeV]}$" )
plt.ylabel( r"$\Delta_{l,s}$" )
plt.xlim( 125 , 220 )
plt.ylim( 0 , 0.95 )

def Draw( ax , x , y , ye , color , ls , label ):
    ax.fill_between( x , y-ye , y+ye , alpha=0.2 , color=color )
    ax.plot( x , y , color=color , ls=ls , label=label )
    return

Draw( plt  , dat[:,0] , dat[:,7] , dat[:,8]*.01 , "blue" , "-" , r"$\Delta_{l,s}$" )

# plot the renormalized Polyakov loop
# plt2 = plt.twinx()
# plt.tick_params()
# plt2.set_ylabel( r"$L_{\rm ren}$" )
# plt2.set_ylim( 0 , 0.475 )

# plt.plot( [0] , [0] , color="blue" , lw=2.5 , ls="-" , label= r"$\Delta_{l,s}$" ) #for label
# Draw( plt2 , dat[:,0] , dat[:,3] , dat[:,4]*.001 ,
#       "red" , "--" , r"$L_{\rm ren}$" )
# plt.plot( [0] , [0] , color="blue" , lw=2.5 , ls="solid" ) 
# plt.legend( loc="upper right" , ncol=2 )

plt.savefig("chiral_lattice.pdf")
