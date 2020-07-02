import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

def cell_model(t, z):
    
    k1=0.04 
    k2p=0.04
    k2pp=1
    k3p=1
    k3pp=10
    k4=35
    k5p=0.005
    k5pp=0.2
    k6=0.1
    k7=1
    k8=0.5
    k9=0.1
    k10=0.02 
    J3=0.04
    J4=0.04
    J5=0.3
    J7=0.001
    J8=0.001
    n=4
    Mad=1
    mu=0.01
    mstar=10

    cycB, cdh1, cdc20t, cdc20a, IEP, m = z

    # d/dt cycB
    eq1 = k1 - (k2p + k2pp * cdh1) * cycB

    # d/dt Cdh1
    eq2 = ((k3p+k3pp*cdc20a)*(1-cdh1))/(J3+1-cdh1)-(k4*m*cycB*cdh1)/(J4+cdh1)

    # d/dt Cdc20t
    eq3 = k5p + k5pp*((cycB*m/J5)**n)/((1 + (cycB*m/J5))**n) - k6*cdc20t 
    
    # d/dt Cdc20a
    eq4 = (k7*IEP*(cdc20t - cdc20a)/(J7 + cdc20t - cdc20a)) - (k8*Mad*cdc20a)/(J8+cdc20a) - k6*cdc20a

    # d/dt IEP
    eq5 = k9*m*cycB*(1 - IEP) - k10*IEP

    # d/dt m
    eq6 = mu*m*(1-m/mstar)

    return [eq1, eq2, eq3, eq4, eq5, eq6]



sol = integrate.solve_ivp(cell_model, [0, 150], [0.6, 0, 1.5, 0.6, 0.6, 1.1],
                          dense_output = True)

t = np.linspace(0, 150, 1500)
z = sol.sol(t)
plt.plot(t, z.T)
plt.legend(['cycB', 'cdh1', 'cdc20t', 'cdc20a', 'IEP', 'm'])
plt.show()
