"""
This module includes different optical beam shapes
"""
import numpy as np
from scipy.special import assoc_laguerre
import my_functions.functions_general as fg


def LG_simple(x, y, z=0, l=1, p=0, width=1, k0=1, x0=0, y0=0):
    """
    Classic LG beam
    :param l: azimuthal index
    :param p: radial index
    :param width: beam waste
    :param k0: wave number
    :param x0: center of the beam in x
    :param y0: center of the beam in y
    :return: complex field
    """

    def laguerre_polynomial(x, l, p):
        return assoc_laguerre(x, p, l)

    x = x - x0
    y = y - y0
    zR = k0 * width ** 2
    E = (np.sqrt(np.math.factorial(p) / (np.pi * np.math.factorial(np.abs(l) + p)))
         * fg.rho(x, y) ** np.abs(l) * np.exp(1j * l * fg.phi(x, y))
         / (width ** (np.abs(l) + 1) * (1 + 1j * z / zR) ** (np.abs(l) + 1))
         * ((1 - 1j * z / zR) / (1 + 1j * z / zR)) ** p
         * np.exp(-fg.rho(x, y) ** 2 / (2 * width ** 2 * (1 + 1j * z / zR)))
         * laguerre_polynomial(fg.rho(x, y) ** 2 / (width ** 2 * (1 + z ** 2 / zR ** 2)), np.abs(l), p)
         )
    return E


def trefoil(x, y, z, w, width=1, k0=1, z0=0., aCoeff=None, coeffPrint=False):
    z = z - z0
    H = 1.0
    if aCoeff is not None:
        aSumSqr = 0.1 * np.sqrt(sum([a ** 2 for a in aCoeff]))
        aCoeff /= aSumSqr
    else:
        a00 = 1 * (H ** 6 - H ** 4 * w ** 2 - 2 * H ** 2 * w ** 4 + 6 * w ** 6) / H ** 6
        a01 = (w ** 2 * (1 * H ** 4 + 4 * w ** 2 * H ** 2 - 18 * w ** 4)) / H ** 6
        a02 = (- 2 * w ** 4 * (H ** 2 - 9 * w ** 2)) / H ** 6
        a03 = (-6 * w ** 6) / H ** 6
        a30 = (-8 * np.sqrt(6) * w ** 3) / H ** 3
        aCoeff = [a00, a01, a02, a03, a30]
        aSumSqr = 0.1 * np.sqrt(sum([a ** 2 for a in aCoeff]))
        aCoeff /= aSumSqr
    if coeffPrint:
        print(aCoeff)
        print(f'a00 -> a01 -> a02 ->... -> a0n -> an0:')
        for i, a in enumerate(aCoeff):
            print(f'a{i}: {a:.3f}', end=',\t')
    field = (aCoeff[0] * LG_simple(x, y, z, l=0, p=0, width=width, k0=k0) +
             aCoeff[1] * LG_simple(x, y, z, l=0, p=1, width=width, k0=k0) +
             aCoeff[2] * LG_simple(x, y, z, l=0, p=2, width=width, k0=k0) +
             aCoeff[3] * LG_simple(x, y, z, l=0, p=3, width=width, k0=k0) +
             aCoeff[4] * LG_simple(x, y, z, l=3, p=0, width=width, k0=k0)
             )
    return field
