# discrete dipole approximation

import numpy as np 
import datetime
from scattermeta.consts import C_CONST, Mu0_CONST, Epsilon0_CONST

class Point:
    def __init__(self, position, radius, epsilon):
        self.position = np.array(position)
        self.radius = radius
        self.epsilon = epsilon

class PointList:
    def __init__(self, verbose=False):
        self.pointlist = list()
        self.verbose = verbose

    def add_point(self, position, radius, epsilon):
        point = Point(position,radius,epsilon)
        self.pointlist.append(point)

    def get_list(self):
        return self.pointlist

class DDA:
    def __init__(self, pointlist, verbose=False):
        self.pointlist = pointlist.get_list()
        self.verbose = verbose

        self.frequency = None
        self.__dipoles = None

    def get_far_field(self, frequency):
        # define at the spherical coordination

        # TODO

        pass

    def get_spherical_field(self, frequency, R=None, Npoints=[100, 100]):
        N = len(self.pointlist)
        positions = np.array([point.postions for point in self.pointlist])
        center = np.average(positions, 0)

        if R is None:
            distances = list(map(np.linalg.norm, positions - center))
            R = 2*np.max(distances)     

        Ntheta, Nphi = Npoints
        phis = np.linspace(0, 2*np.pi, Nphi)
        thetas = np.linspace(0, np.pi, Ntheta)

        electric_fields, angle_list = list(), list()

        # xs, ys, zs = list(), list(), list()
        for i in range(Ntheta):
            for j in range(Nphi):
                x = R*np.sin(thetas[i])*np.cos(phis[j])
                y = R*np.sin(thetas[i])*np.sin(phis[j])
                z = R*np.cos(thetas[i])
                position = np.array([x, y, z])
                angle_list.append(thetas[i], phis[j])

                electric_field = self.get_electric_field(position, frequency)
                electric_fields.append(electric_field)

        #         xs.append(x)
        #         ys.append(y)
        #         zs.append(z)

        # xs = np.array(xs)
        # ys = np.array(ys)
        # zs = np.array(zs)

        electric_fields = np.array(electric_field)
        angle_list = np.array(angle_list)
        return angle_list, electric_fields

    def get_electric_field(self, position):  
        if self.__dipoles is None:
            raise Exception("Null solution! Please run the function calculate firstly!")

        frequency = self.frequency
        position = np.array(position)
        electric_field = np.zeros(3, dtype=complex)
        k0 = 2*np.pi*frequency/C_CONST
        N = len(self.pointlist)
        for i in range(N):
            r1 = position
            r2 = self.pointlist[i].position
            green_function = k0**2/Epsilon0_CONST*self.__green_function(r1, r2, frequency)
            electric_field += green_function*self.__dipoles[3*i:(3*i+1)]

        return electric_field

    def __green_function(self, r1, r2, frequency):
        r1, r2 = np.array(r1), np.array(r2)
        R = np.linalg.norm(r1-r2)
        k0 = 2*np.pi*frequency/C_CONST
        phase_coefficient = np.exp(1j*k0*R)
        
        # TODO Just for the test, need to rewrite the Green function here!

        green_function = Mu0_CONST/(4*np.pi*R)*phase_coefficient

        return green_function

    def __get_polarization(self, epsilon, radius):
        polarization = 4*np.pi*radius**3*Epsilon0_CONST*(epsilon - 1)/(epsilon + 2)
        return polarization

    def __dipole_eqaution(self, frequency):
        k0 = 2*np.pi*frequency/C_CONST
        N = len(self.pointlist)
        coupling_matrix = np.zeros((3*N, 3*N), dtype=complex)
        polarizations = np.zeros(3*N, dtype=complex)

        if self.verbose:
            print("Initiating the coupling dipole matrix...", end="")

        for i in range(N):
            point = self.pointlist[i]
            polarization = self.__get_polarization(
                epsilon=point.epsilon,
                radius=point.radius
                )
            polarizations[3*i:3*(i+1)] = polarization

            for di in range(3):
                coupling_matrix[3*i+di,3*i+di] = 1

        for i in range(N):
            for j in range(N):
                if i != j:
                    pointi = self.pointlist[i]
                    pointj = self.pointlist[j]

                    green_function = self.__green_function(
                        r1=pointi.position,
                        r2=pointj.position,
                        frequency=frequency
                    )

                    # consider isotropic polarization
                    for di in range(3):
                        for dj in range(3): 
                            coupling = -polarizations[3*j + dj]*k0**2/Epsilon0_CONST*green_function
                            coupling_matrix[3*i + di, 3*j + dj] = coupling
        
        if self.verbose:
            print("Done!")

        return coupling_matrix, polarizations

    def calculate(self, frequency, input_electric_field, input_direction):
        if self.verbose:
            starttime = datetime.datetime.now()
            print("Calcuating at the frequency [{frequency} Hz]...".format(
                    frequency=frequency
                ))

        coupling_matrix, polarizations = self.__dipole_eqaution(frequency)

        N = len(self.pointlist)
        k0 = 2*np.pi*frequency/C_CONST
        b = np.zeros(3*N, dtype=complex)

        input_electric_field = np.array(input_electric_field)
        input_electric_field = input_electric_field/np.linalg.norm(input_electric_field)
        input_direction = np.array(input_direction)
        input_direction = input_direction/np.linalg.norm(input_direction)
        k = k0*input_direction

        for i in range(N):
            for axes in range(3):
                position = self.pointlist[i].position
                phase_coeff = np.exp(1j*k.dot(position))
                b[3*i + axes] = polarizations[3*i + axes]*input_electric_field[axes]*phase_coeff

        if self.verbose:
            print("Calculating dipoles... ", end="")

        # print(coupling_matrix)
        # print(np.linalg.det(coupling_matrix))
        # print(b)
        dipoles = np.linalg.solve(coupling_matrix, b)
        self.frequency = frequency
        self.__dipoles = dipoles

        if self.verbose:
            endtime = datetime.datetime.now()
            second = (endtime - starttime).seconds
            minute = int(np.round(second/60))
            second = second - 60*minute
            print("Done! Spend time: [{minute} min {second} sec]".format(
                minute=minute,
                second=second
            ))
        
        
