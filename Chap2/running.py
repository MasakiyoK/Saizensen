
import numpy as np
import matplotlib.pyplot as plt
from plt_setting import *

plt.xscale("log")
#plt.yscale("log")
plt.xlabel( r"$\mu/\Lambda_{\rm QCD}$" )
plt.ylabel( r"$\alpha(\mu)=g^2(\mu)/4\pi$" )
plt.xlim( 0.8 , 800 )
plt.ylim( 0. , 1.3 )
ytics(0.2)
plt.xticks( [1 , 3 , 10 , 30 , 100 , 300 ] ,
            ["1" , "3" , "10" , "30" , "100" , "300" ] )

muLam = np.linspace( 0. , np.log(1000.) , 50 )
muLam = np.exp(muLam)
alpha = 12.* np.pi / ( 33. - 12. ) / np.log( muLam**2 )
plt.plot( muLam , alpha+.1 )
plt.savefig( "running.pdf" )
