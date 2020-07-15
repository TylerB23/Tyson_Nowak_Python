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

    cycB, cdh1, cdc20t, m = z

    # d/dt cycB
    eq1 = k1 - (k2p + k2pp * cdh1) * cycB

    # d/dt Cdh1
    eq2 = ((k3p+k3pp*cdc20t)*(1-cdh1))/(J3+1-cdh1)-((k4*m*cycB*cdh1)/(J4+cdh1))

    # d/dt Cdc20t
    eq3 = k5p + k5pp*((cycB*m/J5)**n)/((1 + (cycB*m/J5))**n) - k6*cdc20t 
    
    # d/dt m
    eq6 = mu*m*(1-m/mstar)
    #scalar = 2 * np.pi / 75
    #eq6 = -1 * scalar * np.sin(scalar * m)

    return [eq1, eq2, eq3, eq6]
'''
Solving
'''
# steady state between cycB and cdh1:
# G1: 0.039, 0.97 This steady state is lost to bifurcation at m = 0.53
# S-G2-M: 0.9, 0.0045
# Note: cdh1 is in (0,1) per the model.  Operates as a switch
# [cycB, cdh1, cdc20t, m]
g1_ss_ic = [0.04, 0.95, .05, 0.4]
sg2m_ss_ic = [0.9, 0.005, 1, 0.9]
cos_ic = [0.15, 1, 0.2, 0]

def event(x, y):
    return y[0] - 0.1 # when cycB crosses +0.1 from above
event.direction = -1
event.terminal = True

def fig5_mod(ic, t_range=[0,150], plot=False):

    sol0 = integrate.solve_ivp(cell_model, t_range, ic,
                               events=[event], dense_output=True)
    print(sol0.t_events)
    time = sol0.t
    y_vals = sol0.y

    stop = sol0.status
    flag = sol0.status
    while(flag and stop < t_range[1]):
        new_ic = y_vals[:,-1]
        # Split the cell
        y_vals[:, -1][3] /= 2

        sol = integrate.solve_ivp(cell_model, [stop, t_range[1]], new_ic,
                                  events=[event], dense_output = True)
        print(sol.t_events)
        
        time = np.append(time, sol.t)
        y_vals = np.append(y_vals, sol.y, axis=1)
        
        stop = sol.t[-1]
        flag = sol.status

    if plot:
        num_timesteps = time.shape[0]
        t = np.linspace(0, t_range[1], num_timesteps)

        plt.subplot(311)
        plt.plot(t, y_vals[3,:])
        plt.legend(['m'])

        plt.subplot(312)
        plt.plot(t, y_vals[0,:])
        plt.plot(t, y_vals[1,:])
        plt.legend(['cycB', 'cdh1'])

        plt.subplot(313)
        plt.plot(t, y_vals[0,:])
        plt.plot(t, y_vals[2,:])
        plt.legend(['cycB', 'cdc20t'])

        plt.show()

    return time, y_vals
