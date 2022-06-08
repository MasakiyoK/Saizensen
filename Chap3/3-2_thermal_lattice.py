### Saizensen, Fig. 3.2: energy density & pressure, lattice result
### June 2022 / Masakiyo Kitazawa

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from plt_setting import *

# lattice result of thermodynamics:
# BW from Phys. Lett. B 730 (2014) [arXiv:1309.5258]
# https://arxiv.org/abs/1309.5258
# HQ from Phys. Rev. D90 , 094503 [arXiv:1407.6387]
# https://arxiv.org/abs/1407.6387
BW = np.loadtxt( "WB-EoS.dat" )
HQ = np.loadtxt( "thermHotQCD.csv" )

def Draw( x , y , ye1 , ye2 , color , linestyle , label="" ):
    """ Plot lattice data; center line and error band """
    plt.fill_between( x , y-ye1 , y+ye2 , alpha=0.2 , color=color )
    plt.plot( x , y , color=color , ls=linestyle , label=label )
    return

# energy density
Draw( BW[:,0] , BW[:,5] , BW[:,6] , BW[:,6] , "brown" , "solid" , "BW" )
Draw( HQ[:,0] , HQ[:,4] , .01*HQ[:,5] , .01*HQ[:,6] , "red" , "dashed" , "HQ" )

# pressure
Draw( BW[:,0] , 3.*BW[:,3] , 3.*BW[:,4] , 3.*BW[:,4] , "green" , "solid" )
Draw( HQ[:,0] , 3.*HQ[:,1] , .003*HQ[:,2] , .003*HQ[:,3] , "blue" , "dashed" )

# Stefan-Boltzmann limit Eq.(3.9)
SB = pi*pi/30. * ( 16. + 7./8.*36 )
plt.plot( [300,500] , [SB,SB] , ls="--" , color="black" , lw=1.5 )

# HRG model Eq.(3.12)
HRG = np.loadtxt( "HRGpe.dat" )
plt.plot( HRG[HRG[:,0]<211,0] , 3.*HRG[HRG[:,0]<211,1] , ls=":" , color="black" , lw=2 )
plt.plot( HRG[HRG[:,0]<196,0] , HRG[HRG[:,0]<196,2] , ls=":" , color="black" , lw=2 )

### plot
plt.xlabel( r"$T \,\, \mathrm{[MeV]}$" )
plt.ylabel( r"$\varepsilon/T^4 \, , \,\, 3p/T^4$" )

plt.xlim( 110 , 410 )
plt.ylim( 0 , 17 )
text( 220 , 11.8 , r"$\varepsilon/T^4$" , fontsize=25 )
text( 250 , 5.7 , r"$3p/T^4$" , fontsize=25 )
text( 275 , 15.2 , r"$\mathrm{SB}$" , fontsize=25 )

plt.legend( loc="lower right" )
plt.savefig( "thermal_lattice.pdf" )
