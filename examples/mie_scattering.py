import sys
sys.path.append('..')

import matplotlib.pyplot as plt 
import numpy as np 

from scattermeta.mie import MieScatter

if __name__ == '__main__':
    N = 1000
    miescatter = MieScatter(r=1.0, n=1.33)
    ldas = np.linspace(0.1, 3, N)
    C_exts = np.zeros(N)
    for i in range(N):
        C_ext, C_sca, C_abs = miescatter.cal_C(ldas[i])
        C_exts[i] = C_ext
    
    plt.plot(ldas, C_exts)
    plt.xlabel("$\lambda$")
    plt.ylabel("$C_{ext}$")
    plt.show()