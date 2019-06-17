# Mie scattering
import numpy as np
import scipy.special as sp

class MieScatter:
    def __init__(self, r, n, nb=1.0, verbose=False):
        self.r = r
        self.n = n
        self.nb = nb

    def __cal_ab(self, lda):
        # calculate the coefficient an, bn
        x = 2*np.pi/lda*self.r
        m = self.n/self.nb
        z = m*x

        n_max = int(np.round(2 + x + 4*x**(1/3)))
        n_mx = int(np.round(np.max([n_max, np.abs(z)]) + 16))
        n = np.arange(1, n_max + 1)
        n_u = n + 0.5

        s_x  = np.sqrt(0.5*np.pi*x)
        p_x = s_x*sp.jv(n_u, x)
        p_1_x = np.zeros(n_max, dtype='complex')
        p_1_x[0] = np.sin(x)
        p_1_x[1:] = p_x[:(n_max-1)]

        ch_x = -s_x*sp.yv(n_u, x)
        ch_1_x = np.zeros(n_max, dtype='complex')
        ch_1_x[0] = np.cos(x)
        ch_1_x[1:] = ch_x[:(n_max-1)]

        gsx = p_x - 1j*ch_x
        gs1x = p_1_x - 1j*ch_1_x

        dnx = np.zeros(n_mx, dtype='complex')
        for i in range(n_mx-1,0,-1):
            dnx[i-1] = i/z - 1/(dnx[i] + i/z)

        dn = dnx[n]

        da = dn/m + n/x 
        db = m*dn + n/x

        an = (da*p_x - p_1_x)/(da*gsx - gs1x)
        bn = (db*p_x - p_1_x)/(db*gsx - gs1x)

        return an, bn


    def cal_C(self, lda):
        an, bn = self.__cal_ab(lda)
        n_max = len(an)
        cn = 2*np.arange(1,n_max+1) + 1
        x = 2*np.pi/lda*self.r
        C = np.sum(cn*np.real(an + bn))
        C_ext = 2*C/(x*x)
        C = np.sum(cn*(np.abs(an) + np.abs(bn)))
        C_sca = 2*C/(x*x)
        C_abs = C_ext - C_sca
        
        return C_ext, C_sca, C_abs

    def E_sca(self, pos):
        pass