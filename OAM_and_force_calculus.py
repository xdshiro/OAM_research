import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np
import matplotlib.pyplot as plt


def W_energy_3D(EArray, xArray=None, yArray=None):
    """
    total power in Oxy plane
    :param EArray:
    :param xArray:
    :param yArray:
    :return:
    """
    if xArray is None or yArray is None:
        shape = np.shape(EArray)
        xArray = np.arange(shape[0])
        yArray = np.arange(shape[1])
    dx = xArray[1] - xArray[0]
    dy = yArray[1] - yArray[0]
    W = np.real(np.sum(np.conj(EArray) * EArray) * dx * dy)
    return W


def Jz_calc_no_conj_3D(EArray, x0=None, y0=None, z0=None, xArray=None, yArray=None, zArray=None):
    EArray = np.array(EArray)
    Er, Ei = np.real(EArray), np.imag(EArray)
    if xArray is None or yArray is None or zArray is None:
        print("xyz are not defined -> done automatically")
        shape = np.shape(EArray)
        xArray = np.arange(shape[0])
        yArray = np.arange(shape[1])
        zArray = np.arange(shape[2])
    if x0 is None:
        x0 = (xArray[-1] + xArray[0]) / 2
    if y0 is None:
        y0 = (yArray[-1] + yArray[0]) / 2
    if z0 is None:
        return (zArray[-1] + zArray[0]) / 2
    x = np.array(xArray)
    y = np.array(yArray)
    z = np.array(zArray)
    dx = xArray[1] - xArray[0]
    dy = yArray[1] - yArray[0]
    sumJz = 0
    for i in range(1, len(xArray) - 1, 1):
        for j in range(1, len(yArray) - 1, 1):
            dErx = (Er[i + 1, j] - Er[i - 1, j]) / (2 * dx)
            dEry = (Er[i, j + 1] - Er[i, j - 1]) / (2 * dy)
            dEix = (Ei[i + 1, j] - Ei[i - 1, j]) / (2 * dx)
            dEiy = (Ei[i, j + 1] - Ei[i, j - 1]) / (2 * dy)
            # dErx = (Er[i + 1, j] - Er[i, j]) / (dx)
            # dEry = (Er[i, j + 1] - Er[i, j]) / (dy)
            # dEix = (Ei[i + 1, j] - Ei[i, j]) / (dx * 2)
            # dEiy = (Ei[i, j + 1] - Ei[i, j]) / (dy)
            # print(x[i] * Er[i, j] * dEiy, - y[j] * Er[i, j] * dEix, -
            #           x[i] * Ei[i, j] * dEry, + y[j] * Ei[i, j] * dErx)
            sumJz += (x[i] * Er[i, j] * dEiy - y[j] * Er[i, j] * dEix -
                      x[i] * Ei[i, j] * dEry + y[j] * Ei[i, j] * dErx)
    # Total moment
    Jz = (sumJz * dx * dy)
    W = W_energy_3D(EArray)
    print(f'Total OAM charge = {Jz / W}\tW={W}')
    return Jz


if __name__ == '__main__':
    xMinMax = 4
    yMinMax = 4
    zMinMax = 1
    zRes = 60
    xRes = yRes = 60
    xyzMesh = fg.create_mesh_XYZ(xMinMax, yMinMax, zMinMax, xRes, yRes, zRes, zMin=None)
    # modes = [(0, 0), (0, 1), (0, 2), (1, 0), (2, 0), (2, 1), (2, 2)]
    # modes = [(0, 0), (0, 1), (1, 0), (1, 1)]
    # lpArray = [0, 0, 0, 0]
    # diapason = [1, 1, 1, 1]
    modes = [(0, 1), (1, 0)]
    coeff = [0, 1]
    beam = bp.LG_combination(*xyzMesh, coefficients=coeff, modes=modes)
    print(fg.Jz_calc_no_conj(beam[:, :, zRes // 2]))
    pl.plot_2D(np.abs(beam[:, :, zRes // 2]))
    plt.show()
