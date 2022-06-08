### Saizensen, Fig. 3.1: hadron resonance gas model, energy density
### June 2022 / Masakiyo Kitazawa

import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from plt_setting import *

# list of hadrons from Phys. Rev. D104 , 074512 
# https://arxiv.org/abs/2107.10011
list = np.loadtxt( "QM_hadron_list_ext_strange_2020.txt" , usecols=(1,6,7) )
print( list.shape )

def E( m , T , BF ):
    """ Energy density of free particle Eq.(3.3); BF=1:boson / BF=-1:fermion """
    def integrand( p ):
        Ep = sqrt( m*m + p*p )
        return p*p * Ep * ( 1./( np.exp( Ep/T ) - BF ) if Ep/T<35. else np.exp( -Ep/T ) )
    # upper: upper limit of p integral estimated from
    # exp( -sqrt(p^2+m^2)/T ) = exp(-100) * exp( -m/T )
    upper = 100. * T * sqrt( 1. + m/T/50. )
    return 1./(2.*np.pi*np.pi) / T**4 * integrate.quad( integrand , 0. , upper )[0]

def E_sum( T , Mmax=3000. ):
    """ Sum of E up to mass Mmax """
    ans = 0.
    for had in list:
        m , g , BF = had[0] , had[1] , had[2]
        if m<Mmax : 
            ans += g * E( m , T , BF )
    return ans

### plot
plt.xlabel( r"$T\, \mathrm{[MeV]}$" )
plt.ylabel( r"$\varepsilon/T^4$" )
plt.xlim( 0 , 156 )
plt.ylim( 0 , 3.6 )
xtics( 50. )

Tlist = np.linspace( 1. , 160. , 51 )
plt.plot( Tlist , [ E_sum( T) for T in Tlist ] , label=r"$\mathrm{HRG}$" )
plt.plot( Tlist , [ E_sum( T , 1000. ) for T in Tlist ] , ls="-." , label=r"$m<1~\mathrm{GeV}$" )
plt.plot( Tlist , [ E_sum( T ,  500. ) for T in Tlist ] , ls="--" , label=r"$\pi,~K$" )
plt.plot( Tlist , [ E_sum( T ,  140. ) for T in Tlist ] , ls=":" , label=r"$\pi$" )

plt.legend()
plt.savefig( "HRG_e.pdf" )
