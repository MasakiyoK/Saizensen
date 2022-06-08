### Saizensen, Fig. 7.2:
### Masses of sigma and pi modes in 2-flavor NJL model in the chiral limit
### June 2022 / Masakiyo Kitazawa

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
from plt_setting import *

from NJL_funcs2 import *

Tc = Tc_m0()
print( "T_c = " , Tc )

def Mq( T ):
    """ constituent quark mass for a given T at m=0 """
    if T>=Tc : return 0.
    return optimize.brentq( GapEq , 0. , 350. , args=(T,0.) )

def M_S( T ):
    """ mass of scalar meson m=0 Eq.(7.37) """
    if T<=Tc : return 2.*Mq( T )
    equation = lambda om : 1. - 2.* Gs * PiS0_Re( om , 0. , T )
    return optimize.brentq( equation , 1.e-2 , 700. )

def M_PS( T ):
    """ mass of pseudo-scalar meson m=0 """
    if T<=Tc : return 0.
    return M_S( T )

### plot
plt.xlabel( r"$T ~{\rm [MeV]}$" )
plt.ylabel( r"${\rm mass~[MeV]}$" )
plt.xlim( 0. , 312. )
plt.ylim( 0. , 720. )

gridT = np.concatenate([ np.linspace( 1. , Tc-20. , 20 ) ,
                         np.linspace( Tc-19. , Tc+19. , 101 ) ,
                         np.linspace( Tc+20. , 320. , 20 ) ])
plt.plot( gridT , [ M_S( T ) for T in gridT ]  , ls='-' , label=r'$\sigma$')
plt.plot( gridT , [ M_PS( T ) for T in gridT ]  , ls='--' , label=r'$\pi$')

plt.fill_between( gridT , [ 2.*Mq( T ) for T in gridT ] , 720. , color='yellow' , alpha=0.2 )

plt.text( 250 , 650 , r"$m=0$" )

plt.legend( loc='center left' )
plt.savefig( "NJL_meson_mass_m0.pdf" )
