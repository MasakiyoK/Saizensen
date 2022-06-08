### Saizensen, Fig. 7.2:
### Masses of sigma and pi modes in 2-flavor NJL model for m=4 MeV
### June 2022 / Masakiyo Kitazawa

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
from plt_setting import *

from NJL_funcs2 import *

m = 4.

def Mq( T , m ):
    """ constituent quark mass for a given T for finite m """
    return optimize.brentq( GapEq , m , 350. , args=(T,m) )

Tc_m0 = Tc_m0() #critical temperature at m=0
print( "T_c (m=0) = " , Tc_m0 )

def M_S( T , m ):
    """ mass of scalar meson for finite m; Eq.(7.37) """
    M_q = Mq(T,m)
    equation = lambda om : 1. - 2.* Gs * PiS0_Re( om , M_q , T )
    return optimize.brentq( equation , 0. , 700. )

def M_PS( T , m ):
    """ mass of pseudo-scalar meson for finite m; Eq.(7.37) """
    M_q = Mq(T,m)
    equation = lambda om : 1. - 2.* Gs * PiPS0_Re( om , M_q , T )
    return optimize.brentq( equation , 1.e-1 , 700. )

# pion mass at T=0
print( "M_pi (T=0) = " , M_PS(1.e-5,m) )

def M_PS0( T ):
    """ mass of pseudo-scalar meson m=0; Eq.(7.37) """
    if T<=Tc_m0 : return 0.
    equation = lambda om : 1. - 2.* Gs * PiPS0_Re( om , 0. , T )
    return optimize.brentq( equation , 1.e-2 , 700. )

### plot
plt.xlabel( r"$T ~{\rm [MeV]}$" )
plt.ylabel( r"${\rm mass~[MeV]}$" )
plt.xlim( 0. , 310. )
plt.ylim( 0. , 720. )

gridT = np.linspace( 0. , 320. , 100 )
plt.plot( gridT , [ M_S( T,m ) for T in gridT ] , ls='-' , label=r'$\sigma$')
plt.plot( gridT , [ M_PS( T,m ) for T in gridT ] , ls='--' , label=r'$\pi$')

gridTc = np.concatenate([ np.linspace( 1. , Tc_m0 , 2 ) ,
                          np.linspace( Tc_m0+.01 , Tc_m0+20. , 40 ) ,
                          np.linspace( Tc_m0+22. , 320. , 20 ) ])
plt.plot( gridTc , [ M_PS0(T) for T in gridTc ] , lw=2 , ls='-.' , label=r'$\pi ~(m=0)$')

plt.plot( gridT , [ 2.*Mq( T,m ) for T in gridT ]  , lw=2 , ls=':' , label=r'$2M$')
plt.fill_between( gridT , [ 2.*Mq( T,m ) for T in gridT ] , 720. , color='yellow' , alpha=0.2 )

plt.legend( loc='center left' )
plt.text( 230 , 650 , r"$m=4~{\rm MeV}$" )
plt.savefig( "NJL_meson_mass.pdf" )
