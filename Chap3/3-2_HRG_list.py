### Saizensen, Fig. 3.2: data production of HRG energy density & pressure
### make "HRGpe.dat"
### See Eq.(3.12)
### June 2022 / Masakiyo Kitazawa

import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from plt_setting import *

# data file from Phys. Rev. D104 , 074512 [arXiv:2107.10011]
# https://arxiv.org/abs/2107.10011
list = np.loadtxt( "QM_hadron_list_ext_strange_2020.txt" , usecols=(1,6,7) )

def E( m , T , BF ):
    """ Energy density of free particle Eq.(3.3); BF=1:boson / BF=-1:fermion """
    def integrand( p ):
        Ep = sqrt( m*m + p*p )
        return p*p * Ep * ( 1./( np.exp( Ep/T ) - BF ) if Ep/T<35. else np.exp( -Ep/T ) )
    upper = 100. * T * sqrt( 1. + m/T/50. )
    return 1./(2.*np.pi*np.pi) / T**4 * integrate.quad( integrand , 0. , upper )[0]

def P( m , T , BF ): 
    """ Pressure of free particle Eq.(3.4); BF=1:boson / BF=-1:fermion """
    def integrand( p ): 
        Ep = sqrt( m*m + p*p )
        return p*p * ( -BF * np.log( 1. -BF * np.exp( -Ep/T ) ) if Ep/T<35. else np.exp( -Ep/T ) )
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

def E_sum( T , Mmax=3000. ):
    """ Sum of E up to mass Mmax """
    ans = 0.
    for had in list:
        m , g , BF = had[0] , had[1] , had[2]
        if m<Mmax : 
            ans += g * E( m , T , BF )
    return ans

Tlist = np.arange( 5. , 251. , 5. )
Plist = [ P_sum( T) for T in Tlist ]
Elist = [ E_sum( T) for T in Tlist ]

list = np.concatenate( [ Tlist , Plist , Elist ] ).reshape( 3 , -1 ).transpose()
np.savetxt( "HRGpe.dat" , list )
