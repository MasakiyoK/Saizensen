### Saizensen, Fig. 3.1: hadron resonance gas model, pressure
### June 2022 / Masakiyo Kitazawa

import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from plt_setting import *

# list of hadrons from Phys. Rev. D104 , 074512 
# https://arxiv.org/abs/2107.10011
list = np.loadtxt( "QM_hadron_list_ext_strange_2020.txt" , usecols=(1,6,7) )
print( list.shape )

def P( m , T , BF ): 
    """ Pressure of free particle Eq.(3.4); BF=1:boson / BF=-1:fermion """
    def integrand( p ): 
        Ep = sqrt( m*m + p*p )
        return p*p * ( -BF * np.log( 1. -BF * np.exp( -Ep/T ) ) if Ep/T<35. else np.exp( -Ep/T ) )
    # upper: upper limit of p integral estimated from
    # exp( -sqrt(p^2+m^2)/T ) = exp(-100) * exp( -m/T )
    upper = 100. * T * sqrt( 1. + m/T/50. )
    return 1./(2.*np.pi*np.pi) / T**3 * integrate.quad( integrand , 0. , upper )[0]

def P_sum( T , Mmax=3000. ):
    """ Sum of P up to mass Mmax """
    ans = 0.
    for had in list:
        m , g , BF = had[0] , had[1] , had[2]
        if m<Mmax : 
            ans += g * P( m , T , BF )
    return ans

### plot
plt.xlabel( "$T ~\mathrm{[MeV]}$" )
plt.ylabel( "$p/T^4$" )

Tlist = np.linspace( 1. , 160. , 51 )

plt.plot( Tlist , [ P_sum( T) for T in Tlist ] , label=r"$\mathrm{HRG}$" )
plt.plot( Tlist , [ P_sum( T , 1000. ) for T in Tlist ] , ls="-." , label=r"$m<1~\mathrm{GeV}$" )
plt.plot( Tlist , [ P_sum( T ,  500. ) for T in Tlist ] , ls="--" , label=r"$\pi,~K$" )
plt.plot( Tlist , [ P_sum( T ,  140. ) for T in Tlist ] , ls=":" , label=r"$\pi$" )

plt.xlim( 0 , 156 )
plt.ylim( 0 , 0.65 )
xtics( 50. )

plt.legend()
plt.savefig( "HRG_p.pdf" )
