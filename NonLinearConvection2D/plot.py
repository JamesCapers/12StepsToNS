import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def main():
    x, y = np.loadtxt("coords.txt", unpack=True)
    Y, X = np.meshgrid(y, x)

    initial_u = np.loadtxt("initial_u.txt")
    initial_v = np.loadtxt("initial_v.txt")
    final_u = np.loadtxt("final_u.txt")
    final_v = np.loadtxt("final_u.txt")

    surfPlot(X, Y, initial_u)
    # colorMapPlot(X, Y, initial)
    plt.title("Initial Condition U")

    surfPlot(X, Y, initial_v)
    # colorMapPlot(X, Y, initial)
    plt.title("Initial Condition V")

    surfPlot(X, Y, final_u)
    # colorMapPlot(X, Y, final)
    plt.title("Final Solution U")

    surfPlot(X, Y, final_v)
    # colorMapPlot(X, Y, final)
    plt.title("Final Solution V")

    plt.show()

def colorMapPlot(X, Y, Z):
    plt.figure()
    plt.pcolor(X, Y, Z, cmap='bwr')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.colorbar()

def surfPlot(X, Y, Z):
    plt.figure()
    ax = plt.gca(projection='3d')
    surf = ax.plot_surface(X, Y, Z, cmap=cm.bwr)
    plt.xlabel("x")
    plt.ylabel("y")

if __name__ == '__main__':
    main()
