import numpy as np 

class Scatter:
    def __init__(self, shape):
        self.shape = shape
        self.structure = np.ones(shape)

    def get_structure(self):
        return np.copy(self.structure)

    def set_structure(self, structure):
        self.structure = np.copy(structure)