import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import optimize

Nf = 2.
Nc = 3.
Lambda = 650.
Gs = 5.01e-6 
Gc = 3.11e-6

func = lambda x ,a,b: a*b

print ( integrate.quad( func , -1. , 5. , args=(1.,2.) , weight="cauchy" , wvar=0.) )
print ( np.log(5))
