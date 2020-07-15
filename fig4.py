import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from gk_func import gk

def fig2_model(t, z):
    
    cycB, cdc20t = z

    cdh1 = gk(k3p + k3pp*cdc20t, k4*m*cycB, j3, J4)

    #d/dt cycB
    eq1 = k1 - (k2p + k2pp * cdh1) * cycB

    #d/dt cdc20t
    eq3 = k5p + k5pp*(cycB*m/J5)**n)/((1+(cycB*m/J5))**n) - k6*cdc20t 

    return [eq1, eq3]

def plot_nc(m, cycB_range=[0,1.2], cdc20t_range=[0,.5]):

    cdc20t_nc = 
