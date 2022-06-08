### Saizensen, Fig.7.3
### NJL model; meson spectral function at mu=0
### Feb. 2022 / Masakiyo Kitazawa

import numpy as np
import matplotlib.pyplot as plt
from plt_setting import *

from NJL_funcs2 import *

# critical temperature
Tc = Tc_m0()
print( "T_c = " , Tc )

def rho( om , T ):
    """ Spectral function """
    Re = 1. - 2.* Gs * PiS0_Re( om , 0. , T )
    Im = -2.* Gs * 2. * NfNc/np.pi * om/2. * ( 1. - 2. * Fermi( om/2. , T ) )
    return -Im / ( Re*Re + Im*Im ) / np.pi

### plot
plt.xlim( 0. , 182. )
plt.ylim( 0. , 1500. )
plt.xlabel( r"$\omega ~{\rm [MeV]}$" )
plt.ylabel( r"$\rho(\omega)$" )

grid = np.linspace( 1. , 200. , 500 )
plt.plot( grid , [ rho( om,1.02*Tc ) for om in grid ]  , ls='-' , label=r'$T=1.02T_c$')
plt.plot( grid , [ rho( om,1.05*Tc ) for om in grid ]  , ls='--' , label=r'$T=1.05T_c$')
plt.plot( grid , [ rho( om,1.1*Tc ) for om in grid ]  , ls=':' , label=r'$T=1.1T_c$')

plt.legend()
plt.savefig( "NJL_meson_spc.pdf" )
