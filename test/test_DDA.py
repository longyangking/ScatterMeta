import sys
sys.path.append('..')
sys.path.append('.')

import numpy as np
from scattermeta.dda import DDA, PointList

if __name__ == "__main__":
    pointlist = PointList()
    dx = 0.1
    for i in range(3):
        for j in range(3):
            pointlist.add_point([-0.1 + dx*i, -0.1 + dx*j,0], 0.02, 2)

    dda = DDA(pointlist, verbose=True)
    input_E_field = [0, 0, 1]
    input_direction = [1, 0, 0]
    dda.calculate(3e5, input_E_field, input_direction)

    Nx, Ny = 100, 100
    xs = np.linspace(-1, 1, Nx)
    ys = np.linspace(-1, 1, Ny)
    field = np.zeros((Nx, Ny), dtype=complex)

    for i in range(Nx):
        for j in range(Ny):
            position = [xs[i], ys[j], 0]
            electric_field = dda.get_electric_field(position)
            Ex, Ey, Ez = electric_field
            field[i, j] = Ez

    import matplotlib.pyplot as plt
    plt.imshow(np.real(field))
    plt.show()
    