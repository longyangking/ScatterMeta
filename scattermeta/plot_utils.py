import numpy as np 
import matplotlib.pyplot as plt 

def plot1d(x, y, title, xlabel, ylabel, xlim. ylim, filename):
    plt.figure()
    plt.plot(x,y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.savefig(filename)

