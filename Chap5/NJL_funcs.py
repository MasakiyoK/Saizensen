### Saizensen, Fig. 5.3-6:
### 2-flavor NJL model
### thermodynamic potential Eq.(5.69)
### June 2022 / Masakiyo Kitazawa

import numpy as np
from scipy import integrate
from scipy import optimize

Lambda = 631. # Eq.(5.71)
g = 5.5e-6 # Eq.(5.71)
NfNc_pi2 = 6. / np.pi / np.pi

def omega_T0( M , m , mu ):
    """ thermodynamic potential at T=0 """
    integrand = lambda p,mass : p*p * ( np.sqrt( mass*mass + p*p ) - mu )
    ans = integrate.quad( integrand , 0. , Lambda , args=(M) )[0] 
    pF = np.sqrt( mu*mu - M*M ) if mu>abs(M) else 0. # Fermi momentum
    ans -= integrate.quad( integrand , 0. , pF , args=(M) )[0] 

    #subtract \sigma=0
    ans -= integrate.quad( integrand , 0. , Lambda , args=(m) )[0] 
    pF = np.sqrt( mu*mu - m*m ) if mu>m else 0. # Fermi momentum
    ans += integrate.quad( integrand , 0. , pF , args=(m) )[0] 

    return -NfNc_pi2 * ans + ( M-m )*( M-m ) / ( 4.* g )

def omega( M , m , T , mu , eps=1. ):
    """ thermodynamic potential; Eq. (5.69) """
    if( T==0 ) : return omega_T0( M , m , mu )
    def integrand( p ):
        Ep = np.sqrt( M*M + p*p )
        Epm , Epp = Ep - mu , Ep + mu
        ans = Ep + T * np.log( ( 1. + np.exp( -Epm/T ) )*( 1. + np.exp( -Epp/T ) ) )
        #subtract \sigma=0
        Ep = np.sqrt( m*m + p*p )
        Epm , Epp = Ep - mu , Ep + mu
        ans -= Ep + T * np.log( ( 1. + np.exp( -Epm/T ) )*( 1. + np.exp( -Epp/T ) ) )
        return p*p * ans
    integ = integrate.quad( integrand , 0. , Lambda , epsabs=eps )[0]
    return -NfNc_pi2 * integ + ( M-m )*( M-m ) / ( 4.* g )

def M_T0( mu , m ):
    """ global minimum of omega at T=0 """
    min1 = optimize.minimize_scalar( omega_T0 , bracket=(10.,50.) , args=(m,mu) )
    min2 = optimize.minimize_scalar( omega_T0 , bracket=(280.,320.) , args=(m,mu) )
    return min1.x if min1.fun < min2.fun else min2.x
