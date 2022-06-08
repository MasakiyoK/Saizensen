### Saizensen, Fig. 5.9:
### Phase diagram of 2-flavor NJL model with the CSC phase
### Feb. 2022 / Masakiyo Kitazawa

import numpy as np
from scipy import integrate
from scipy import optimize
import matplotlib as plt
from pylab import *
from plt_setting import *

Lambda = 631.
Gs = 5.5e-6
Gd = .7 * Gs
Nf_pi2 = 2. / np.pi / np.pi
Tmin = 1.
Err = 1.e-1

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

def Minimum2D( T , mu ):
    """ Search for global minimum of omega starting from two initial points """
    function = lambda P : omega( P[0] , P[1] , T , mu )
    min1 = optimize.minimize( function , [280.,0.] , method='Nelder-Mead' )
    min2 = optimize.minimize( function , [0.,50.] , method='Nelder-Mead' )
    return min1.x if min1.fun < min2.fun else min2.x

def mu_c( T ):
    """ critical mu where M vanishes """
    M_at_minimum = lambda mu : Minimum2D( T , mu )[0] - 1.e-2
    return optimize.bisect( M_at_minimum , 260. , 330. , xtol=Err )

def Tc_chiral( mu ):
    """ Critical temperature of chiral transition
        Assume \Delta=0 and determine T where M vanishes """
    func = lambda T : optimize.minimize_scalar( omega , args=(0.,T,mu) ).x - 1.e-2
    return optimize.bisect( func , 10. , 200. , xtol=Err )

def Tc_CSC( mu ):
    """ Critical temperature of CSC transition
        Assume M=0 and determine T where \Delta vanishes """
    func = lambda T : optimize.minimize_scalar( lambda Delta : omega( 0. , Delta , T , mu ) ).x - 1.e-2
    return optimize.bisect( func , 10. , 80. , xtol=Err )

# use the value of T_crit obtained in 5-6_NJLphaseDchiral.py
T_crit = 84.2
mu_crit = mu_c( T_crit )


print( "Plotting Phase diagram:" )
print( "Calculating chiral critical T..." )
grid = np.linspace( 0. , mu_crit , 15 )
plt.plot( grid , [ Tc_chiral(mu) for mu in grid ] , ls="--" , label="second order")

print( "Calculating CSC critical T..." )
grid = np.linspace( 318. , 550. , 15 )
plt.plot( grid , [ Tc_CSC(mu) for mu in grid ]  , ls="--" )

print( "Calculating critical mu of 1st transition..." )
grid = np.linspace( Tmin , T_crit , 10 )
plt.plot( [ mu_c(T) for T in grid ] , grid , lw=4 , label="first order")
plt.plot( [mu_c(T_crit)] , [T_crit] , marker="o" , ms=10 , label="tri-critical point" , ls="none" , color="C1" )

plt.xlabel( "$\mu_q ~\mathrm{[MeV]}$" )
plt.ylabel( "$T ~\mathrm{[MeV]}$" )
plt.xlim( 0. , 525. )
plt.ylim( 0. , 205. )

plt.annotate( "$M\\ne0$\n$\Delta=0$" ,xy = (140, 50) , fontsize=23 )
plt.annotate( "$M=0$\n$\Delta=0$" ,xy = (350, 85) , fontsize=23 )
plt.annotate( "$M=0$\n$\Delta\\ne0$" ,xy = (390, 15) , fontsize=23 )

plt.legend()
plt.savefig("phaseD_CSC.pdf")
