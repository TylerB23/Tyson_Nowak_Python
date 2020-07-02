import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

def fig4(cdc20t_range=[0,1], accuracy=1000):

    # cycB in terms of cdh1
    cycB_nc = lambda cdh1 : 0.04/(1 + cdh1)


