### Saizensen, Sec. 7
### NJL model; gap equation, meson polarization func
### at mu=0
### June 2022 / Masakiyo Kitazawa

import numpy as np
from scipy import integrate
from scipy import optimize

Gs = 5.5e-6
Lambda = 631.
NfNc = 6.
coef = 2. * NfNc / ( 2. * np.pi*np.pi )

def Fermi( E , T ):
    """ Fermi distribution function """
    if T==0 : return 0. if E>0 else 1.
    if E/T > 35. : return np.exp( -E/T )
    if E/T <-35. : return 1.
    return 1./ ( np.exp( E/T ) + 1. )

def GapEq( M , T , m ):
    """ LHS of gap Equation Eq.(5.70) at \mu=0 """
    def integrand( p ):
        Ep = np.sqrt( p**2 + M**2 )
        return p*p / Ep * ( 1. - 2.* Fermi( Ep,T ) )
    integ = integrate.quad( integrand , 0. , Lambda  )[0]
    return -coef * integ + ( ( M - m )/M if M>0. else 1. ) / ( 2.*Gs )

def PiS0_Re( om , M , T ):
    """ Scalar polarization function Re \Pi_S^R( p=0 , om ); Eq.(7.22) """
    def integrand1( p ):
        Ep = np.sqrt( p**2 + M**2 )
        coef = p * p * p  * ( 1. if M==0 else p / Ep )
        return coef / ( Ep*Ep - om*om/4. ) * ( 1. - 2.* Fermi( Ep,T ) )
    if abs(om) < 2.* M :
        return coef * integrate.quad( integrand1 , 0. , Lambda )[0]
        
    p_pole = np.sqrt( om*om/4. - M*M )
    def integrand2( p ):
        return integrand1( p ) * ( p - p_pole )
    return coef * integrate.quad( integrand2 , 0. , Lambda , weight="cauchy" , wvar=p_pole )[0]

def PiPS0_Re( om , M , T ):
    """ Scalar polarization function Re \Pi_S^R( p=0 , om ); Eq.(7.24) """
    def integrand1( p ):
        Ep = np.sqrt( p**2 + M**2 )
        return p * p * Ep / ( Ep*Ep - om*om/4. ) * ( 1. - 2.* Fermi( Ep,T ) )
    if abs(om) < 2.* M :
        return coef * integrate.quad( integrand1 , 0. , Lambda )[0]
        
    p_pole = np.sqrt( om*om/4. - M*M )
    def integrand2( p ):
        return integrand1( p ) * ( p - p_pole )
    return coef * integrate.quad( integrand2 , 0. , Lambda , weight="cauchy" , wvar=p_pole )[0]

def Tc_m0():
    """ critical temperature at m=0 from Thouless criterion; see page 123 """
    return optimize.brentq( lambda T : GapEq(0.,T,0.) , 100. , 250. )
