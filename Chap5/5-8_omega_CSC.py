### Saizensen, Fig. 5.8:
### 2-flavor NJL model
### 2D contour plot of thermodynamic potential 
### June 2022 / Masakiyo Kitazawa

import numpy as np
from scipy import integrate
from scipy import optimize
import matplotlib.pyplot as plt
from plt_setting import *

Lambda = 631.
Gs = 5.5e-6
Gd = .7 * Gs
Nf_pi2 = 2. / np.pi / np.pi

def omega( M , Delta , T , mu ):
    """ Thermodynamic potential as a function of (M,Delta) Eq.(5.77) """
    def integrand( p ):
        Ep = np.sqrt( M*M + p*p )
        xip = Ep + mu
        xim = Ep - mu
        epsp = np.sqrt( xip*xip + Delta*Delta )
        epsm = np.sqrt( xim*xim + Delta*Delta )
        if( xim<0. ): epsm *= -1.

        arg = ( 1.+np.exp( -xim/T ) ) * ( 1.+np.exp( -xip/T ) )
        tmp = Ep + T * np.log( arg )
        arg = ( 1.+np.exp( -epsm/T ) ) * ( 1.+np.exp( -epsp/T ) )
        tmp += epsm + epsp + 2.* T * np.log( arg )
        return p*p * tmp
    integ = integrate.quad( integrand , 0. , Lambda )[0]
    return M*M / 4./Gs + Delta*Delta / 4./Gd - Nf_pi2 * integ

def Minimum( T , mu ):
    """ Search for local minimum of omega starting from (M,D)=(250.,0.) """
    function = lambda P : omega( P[0] , P[1] , T , mu )
    return optimize.minimize( function , [250.,0.] , method='Nelder-Mead' ).x

### plot
x_max = 365.
y_max = 105.

def Plot( T , mu ):
    print( "T=" , T , " / mu_q=" , mu )
    plt.clf()
    plt.xlabel( r"$M ~\mathrm{[MeV]}$" )
    plt.ylabel( r"$\Delta ~\mathrm{[MeV]}$" )
    plt.title( r"$T=" + str(T) + "~\mathrm{[MeV]} ,~~ \mu_q=" + str(mu) + "~\mathrm{[MeV]}$" )
    plt.xlim( 0. , x_max )
    plt.ylim( 0. , y_max )

    # contour plot
    Mgrid = np.linspace( 0. , x_max , 55 )
    Dgrid = np.linspace( 0. , y_max , 55 )
    Omesh = [[ omega(M,Delta,T,mu) for M in Mgrid ] for Delta in Dgrid ]
    plt.contour( Mgrid , Dgrid , Omesh ,  40 , linewidths=1.5 )

    # mark at the minimum
    Mmin , Dmin = Minimum(T,mu)
    print( "minimum = " , Mmin , Dmin )
    plt.plot( [Mmin] , [Dmin] , marker="x" , ms=15 , mew=5 , clip_on=False)

    plt.savefig( "omegaCSC_T" + str(round(T)) + "mu" + str(round(mu)) + ".pdf" )
    
Plot( 5.  , 350. )
Plot( 50. , 350. )
Plot( 5.  , 300. )
