import numpy as np
import scipy.sparse as sparse

class FDFD:
    def __init__(self, structure, size_unit, bk_material=1.0, verbose=False):
        self.structure = np.array(structure)
        self.shape = self.structure.shape

        self.size_unit = size_unit
        self.bk_material = bk_material
        self.verbose = verbose

    def set_size_unit(self, size_unit):
        self.size_unit = size_unit

    def get_size_unit(self):
        return self.size_unit

    def get_bk_material(self):
        return self.bk_material

    def set_bk_material(self, bk_material):
        self.bk_material = bk_material

    def set_structure(self, structure):
        self.structure = np.copy(structure)

    def get_structure(self):
        return np.copy(self.structure)

    def __get_matrix(self):
        dx = self.dx
        dy = self.dy
        c = 1.0
        epsilon = self.epsilon
        Nx = self.Nx
        Ny = self.Ny
        M  = Nx*Ny
        
        invEPS = sparse.spdiags(c**2/np.reshape(epsilon,M),0,M,M,format='lil')

        DEX = sparse.lil_matrix((M, M), dtype=complex)
        DEX = DEX + sparse.diags([-1],[0],shape=(M,M),format='lil',dtype=complex)
        DEX = DEX + sparse.diags([1],[1],shape=(M,M),format='lil',dtype=complex)

        for ny in range(Ny-1):
            neq = Nx*ny + Nx - 1
            DEX[neq,neq+1] = 0

        dpx = np.exp(-1j*kx)
        for ny in range(Ny):
            neq = Nx*ny + Nx-1
            nv = Nx*ny
            DEX[neq,nv] = dpx
        DEX = DEX/dx

        DEY = sparse.lil_matrix((M, M), dtype=complex)
        DEY = DEY + sparse.diags([-1],[0],shape=(M,M),format='lil',dtype=complex)
        DEY = DEY + sparse.diags([1],[Nx],shape=(M,M),format='lil',dtype=complex)
        dpy = np.exp(-1j*ky)
        for nx in range(Nx):
            neq = Nx*(Ny-1) + nx
            nv = nx
            DEY[neq,nv] = dpy
        DEY = DEY/dy

        DHX = sparse.lil_matrix((M, M), dtype=complex)
        DHX = DHX + sparse.diags([1],[0],shape=(M,M),format='lil',dtype=complex)
        DHX = DHX + sparse.diags([-1],[-1],shape=(M,M),format='lil',dtype=complex)

        for ny in range(1,Ny):
            neq = Nx*ny
            DHX[neq,neq-1] = 0

        dpx = np.exp(1j*kx)
        for ny in range(Ny):
            neq = Nx*ny
            nv = Nx*ny + Nx - 1
            DHX[neq,nv] = -dpx
        DHX = DHX/dx

        DHY = sparse.lil_matrix((M, M), dtype=complex)
        DHY = DHY + sparse.diags([1],[0],shape=(M,M),format='lil',dtype=complex)
        DHY = DHY + sparse.diags([-1],[-Nx],shape=(M,M),format='lil',dtype=complex)
        dpy = np.exp(1j*ky)
        for nx in range(Nx):
            neq = nx
            nv = Nx*(Ny-1) + nx
            DHY[neq,nv] = -dpy
        DHY = DHY/dy

        A = -(DHX.dot(DEX) + DHY.dot(DEY))
        Hamiltonian = invEPS.dot(A)
        return Hamiltonian

    def start(self, angle=0):
        