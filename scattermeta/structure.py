import numpy as np 

class PointCloud:
    def __init__(self, ds, Ns):
        self.ds = ds
        self.Ns = Ns  

        self.cloud = np.zeros(Ns)

    def get_cloud(self):
        return self.mesh

    def get_ds(self):
        return self.ds

    def reset(self, ds, Ns):
        

class Structure:
    def __init__(self, mesh, verbose=False):
        self.shape = shape
        self.dielectric_list = np.ones(shape)

        self.verbose = verbose

    def add_block(self):
        return np.copy(self.structure)
    
    def add_sphere(self):

    def add_cylinder(self):

    def add_ellipsolid(self):

    def get_mesh(self):

    def get_dielectric_list(self):

    def set_dielectric_list(self):
