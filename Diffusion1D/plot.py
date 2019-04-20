import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def main():
    data_files = getDataFiles(sys.argv)
    print("Data files found: ", data_files)
    for f in data_files:
        plotFileXY(f)

    plt.title("Diffusion")
    plt.show()

# Extract list of date files from command line arguments
def getDataFiles(arguments):
    file_list = []
    k = 1 # Ignore script name
    while k < len(arguments):
        filename, file_extension = os.path.splitext(arguments[k])
        if file_extension == ".txt":
            file_list.append(arguments[k])
        k += 1
    return file_list

def plotFileXY(fname):
    xv, yv = np.loadtxt(fname, unpack=True)
    plt.plot(xv, yv, label=fname)
    plt.xlabel("x")
    plt.ylabel("u")
    plt.legend()

if __name__ == '__main__':
    main()
