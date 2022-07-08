"""
this module includes all the general functions used in another modules
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CloughTocher2DInterpolator
from scipy import integrate
import scipy.io as sio


def rho(*r):
    """
    return abs of the vector r
    :param r: [x1, x2, x3...]
    :return: |r|
    """
    ans = 0
    for x in r:
        ans += x ** 2
    return np.sqrt(ans)


def phi(x, y):
    """
    angle phi in the plane
    """
    return np.angle(x + 1j * y)


def distance_between_points(point1, point2):
    """
    distance between 2 points in any dimensions
    :param point1: [x1, ...]
    :param point2: [x2, ...]
    :return: geometrical distance
    """
    deltas = np.array(point1) - np.array(point2)
    return rho(deltas)



def create_mesh_XYZ(xMax, yMax, zMax, xRes=50, yRes=50, zRes=50,
                    xMin=None, yMin=None, zMin=None, indexing='ij', **kwargs):
    """
    creating the mesh using np.meshgrid
    :param xMax: [xMin, xMax] are the boundaries for the meshgrid along x
    :param yMax: [yMin, yMax] are the boundaries for the meshgrid along y
    :param zMax: [zMin, zMax] are the boundaries for the meshgrid along z
    :param xRes: resolution along x
    :param yRes: resolution along y
    :param zRes: resolution along z
    :param xMin: xMin=xMax by default.
    :param yMin: yMin=yMax by default.
    :param zMin: zMin=zMax by default.
    :param indexing: ij is the classic matrix (0,0) left top
    :return: mesh
    """
    if xMin is None:
        xMin = -xMax
    if yMin is None:
        yMin = -yMax
    if zMin is None:
        zMin = -zMax
    xArray = np.linspace(xMin, xMax, xRes)
    yArray = np.linspace(yMin, yMax, yRes)
    zArray = np.linspace(zMin, zMax, zRes)
    return np.meshgrid(xArray, yArray, zArray, indexing=indexing, **kwargs)


def create_mesh_XY(xMax, yMax, xRes=50, yRes=50,
                   xMin=None, yMin=None, indexing='ij', **kwargs):
    """
    creating the mesh using np.meshgrid
    :param xMax: [xMin, xMax] are the boundaries for the meshgrid along x
    :param yMax: [yMin, yMax] are the boundaries for the meshgrid along y
    :param xRes: resolution along x
    :param yRes: resolution along y
    :param xMin: xMin=xMax by default.
    :param yMin: yMin=yMax by default.
    :param indexing: ij is the classic matrix (0,0) left top
    :return: mesh
    """
    if xMin is None:
        xMin = -xMax
    if yMin is None:
        yMin = -yMax
    xArray = np.linspace(xMin, xMax, xRes)
    yArray = np.linspace(yMin, yMax, yRes)
    return np.meshgrid(xArray, yArray, indexing=indexing, **kwargs)


def interpolation_real(field, xArray=None, yArray=None, **kwargs):
    """
    Interpolation of any real 2d matrix into the function
    :param field: initial Real 2D array
    :param xArray: x interval (range)
    :param yArray: y interval (range)
    :param kwargs: extra parameters for CloughTocher2DInterpolator
    :return: CloughTocher2DInterpolator return
    """
    xResolution, yResolution = np.shape(field)
    if xArray is None:
        xArray = list(range(xResolution))
    if yArray is None:
        yArray = list(range(yResolution))
    xArrayFull = np.zeros(xResolution * yResolution)
    yArrayFull = np.zeros(xResolution * yResolution)
    fArray1D = np.zeros(xResolution * yResolution)
    for i in range(xResolution * yResolution):
        xArrayFull[i] = xArray[i // yResolution]
        yArrayFull[i] = yArray[i % yResolution]
        fArray1D[i] = field[i // yResolution, i % yResolution]
    return CloughTocher2DInterpolator(list(zip(xArrayFull, yArrayFull)), fArray1D, **kwargs)


# function interpolate complex 2D array of any data into the function(x, y)
def interpolation_complex(field, xArray=None, yArray=None, fill_value=False):
    """
    function interpolate complex 2D array of any data into the function(x, y)
    :param field: initial complex 2D array
    :param xArray: x interval (range)
    :param yArray: y interval (range)
    :return: Real CloughTocher2DInterpolator, Imag CloughTocher2DInterpolator
    """
    fieldReal = np.real(field)
    fieldImag = np.imag(field)
    real = interpolation_real(fieldReal, xArray, yArray, fill_value=fill_value)
    imag = interpolation_real(fieldImag, xArray, yArray, fill_value=fill_value)

    def f(x, y):
        return real(x, y) + 1j * imag(x, y)

    return f
    # return interpolation_real(fieldReal, xArray, yArray), interpolation_real(fieldImag, xArray, yArray)


def integral_of_function_1D(integrandFunc, x1, x2, epsabs=1.e-5, maxp1=50, limit=50, **kwargs):
    """
    scipy.integrate can only work with real numbers so this function splits the integrand to imaginary and real
    parts and integrates the separately, then combine the answers together
    :param integrandFunc: integrand
    :param x1: lower limit
    :param x2: upper limit
    :return: integral value, (real error, imag error)
    """

    def real_f(x):
        return np.real(integrandFunc(x))

    def imag_f(x):
        return np.imag(integrandFunc(x))

    real_integral = integrate.quad(real_f, x1, x2, epsabs=epsabs, maxp1=maxp1, limit=limit, **kwargs)
    imag_integral = integrate.quad(imag_f, x1, x2, epsabs=epsabs, maxp1=maxp1, limit=limit, **kwargs)
    return real_integral[0] + 1j * imag_integral[0], (real_integral[1:], imag_integral[1:])


def arrays_from_mesh(mesh):
    """
    Functions returns the tuple of x1Array, x2Array... of the mesh
    :param mesh: no-sparse mesh, for 3D: [3][Nx, Ny, Nz]
    :return: for 3D: xArray, yArray, zArray
    """
    xList = []
    for i, m in enumerate(mesh):
        row = [0] * len(np.shape(m))
        row[i] = slice(None, None)
        xList.append(m[tuple(row)])
    xTuple = tuple(xList)
    return xTuple


def reading_file_mat(fileName, fieldToRead="p_charges", printV=False):
    """
    function read the mat file and conver 1 of its fields into numpy array
    :param fileName:
    :param fieldToRead: which field to conver (require the name)
    :param printV: if you don't know the name, set it True
    """
    matFile = sio.loadmat(fileName, appendmat=False)
    if fieldToRead not in matFile:
        printV = True
    if printV:
        print(matFile)
        exit()
    return np.array(matFile[fieldToRead])


def random_list(values, diapason, diapason_complex=None):
    """
    Function returns values + random * diapason
    :param values: we are changing this values
    :param diapason: to a random value from [value - diap, value + diap]
    :return: new modified values
    """
    import random
    if diapason_complex is not None:
        answer = [x + random.uniform(-d, +d) + 1j * random.uniform(-dC, +dC) for x, d, dC
                  in zip(values, diapason, diapason_complex)]
    else:
        answer = [x + random.uniform(-d, +d) for x, d in zip(values, diapason)]
    return answer


if __name__ == '__main__':
    import my_functions.beams_and_pulses as bp
    import my_functions.plotings as pl

    xMinMax, yMinMax = 3, 3
    xRes = yRes = 50
    xyMesh = create_mesh_XY(xMinMax, yMinMax, xRes, yRes)
    beam = bp.LG_simple(*xyMesh, l=2) + bp.LG_simple(*xyMesh, l=1)
    # fieldFunction[0] + 1j * fieldFunction[1]
    xArray = np.linspace(-xMinMax, xMinMax, xRes)
    yArray = np.linspace(-yMinMax, yMinMax, yRes)
    beamInterpolated = interpolation_complex(beam, xArray, yArray)
