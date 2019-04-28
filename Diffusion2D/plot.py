import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def main():
    x, y = np.loadtxt("coords.txt", unpack=True)
    Y, X = np.meshgrid(y, x)

    initial = np.loadtxt("initial.txt")
    final = np.loadtxt("final.txt")

    surfPlot(X, Y, initial)
    # colorMapPlot(X, Y, initial)
    plt.title("Initial Condition")

    surfPlot(X, Y, final)
    # colorMapPlot(X, Y, final)
    plt.title("Final Solution")

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
