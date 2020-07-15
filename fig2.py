import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from gk_func import gk

def fig2(t, z, k1, k2p, k2pp, k3p, k3pp, k4, A, j3, j4, m):
    cycB, cdh1 = z

    eq1 = k1 - (k2p + k2pp * cdh1) * cycB
    eq2 = ((k3p + k3pp * A) * (1-cdh1))/(j3 + 1 - cdh1) - (k4*m*cycB*cdh1)/(j4+cdh1)
    return [eq1, eq2]

def plot_sol(tmin, tmax, cycB0, cdh10):
    args = (0.04, 0.04, 1, 1, 10, 35, 0, 0.04, 0.04, 0.3)
    sol = integrate.solve_ivp(fig2, [0, 100], [0.1, 0.36],
                              args=args, dense_output = True)

    t = np.linspace(0, 100, 1000)
    z = sol.sol(t)
    plt.plot(t, z.T)
    plt.legend(['cycB', 'cdh1'])
    plt.show()

def plot_phase(m=0.3, cycB_range=[0,10], cdh1_range=[0,1], num_grid_points=10):
    # Plot the vector field
    eq1 = lambda cycB, cdh1 : 0.04 - (0.04 + 1 * cdh1) * cycB
    eq2 = lambda cycB, cdh1 : 0 - (35 * m *cycB*cdh1)/(0.04+cdh1)

    cycB_ = np.linspace(cycB_range[0], cycB_range[1], num_grid_points)
    cdh1_ = np.linspace(cdh1_range[0], cdh1_range[1], num_grid_points)

    grid = np.meshgrid(cycB_, cdh1_)
    dfmat = np.zeros((num_grid_points, num_grid_points, 2))
    for nx in range(num_grid_points):
        for ny in range(num_grid_points):
            ddt_cycB = eq1(nx, ny)
            ddt_cdh1 = eq2(nx, ny)
            dfmat[nx, ny, 0] = ddt_cycB
            dfmat[nx, ny, 1] = ddt_cdh1

    plt.quiver(grid[0], grid[1], dfmat[:, :, 0], dfmat[:, :, 1])
    #CS = plt.contour(grid[0], grid[1], dfmat[:, :, 0], colors='r')
    #CS1 = plt.contour(grid[0], grid[1], dfmat[:, :, 1], colors='b')
    #plt.clabel(CS)
    #plt.clabel(CS1)
    plt.show()

    return
def plot_nc(cycB_range=[0.01,1], accuracy=100, m = 0.3): 
    # Plot nullclines
    cycB_nc = lambda cycB : (0.04-cycB*0.04)/(cycB*1)

    cycB_vals = np.linspace(cycB_range[0], cycB_range[1], accuracy)
    plt.plot(cycB_vals, cycB_nc(cycB_vals))

    p = 1/(35*m)
    J3 = 0.04
    J4 = 0.04
    plt.plot(cycB_vals, gk(p, cycB_vals, J3, J4))

    # plot with m = 0.3
    m = 0.3
    p = 1/(35*m)
    plt.plot(cycB_vals, gk(p, cycB_vals, J3, J4))


    # Next, plot with m = 0.6
    m = 0.6
    p = 1/(35*m)
    plt.plot(cycB_vals, gk(p, cycB_vals, J3, J4))

    plt.ylim(0,1)
    plt.xscale('log')
    plt.title('Figure 2 Reproduction')
    plt.xlabel('cycB (log scale)')
    plt.ylabel('cdh1')
    plt.legend(['cycB', 'cdh1, m=user input',' cdh1, m=0.3', 'cdh1, m=0.6'])
    plt.show()

    return
