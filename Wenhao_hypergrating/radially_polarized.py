import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np
from scipy.special import assoc_laguerre

def LG_radially(x, y, z=0, l=1, p=0, width=1, k0=1, x0=0, y0=0, z0=0):
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
    z = z - z0
    zR = k0 * width ** 2
    E = (np.sqrt(np.math.factorial(p) / (np.pi * np.math.factorial(np.abs(l) + p)))
         * fg.rho(x, y) ** np.abs(l) * np.exp(1j * l * fg.phi(x, y))
         / (width ** (np.abs(l) + 1) * (1 + 1j * z / zR) ** (np.abs(l) + 1))
         * ((1 - 1j * z / zR) / (1 + 1j * z / zR)) ** p
         * np.exp(-fg.rho(x, y) ** 2 / (2 * width ** 2 * (1 + 1j * z / zR)))
         * laguerre_polynomial(fg.rho(x, y) ** 2 / (width ** 2 * (1 + z ** 2 / zR ** 2)), np.abs(l), p)
         ) * np.cos(fg.phi(x, y))
    return E


if __name__ == '__main__':
    plot_trefoil = True
    if plot_trefoil:
        xyMesh = fg.create_mesh_XY_old(4.5, 4.5, 150, 150)
        beam = LG_radially(*xyMesh)
        beam = fg.propagator_split_step_3D(beam, dz=1, xArray=None, yArray=None, zSteps=350, n0=1, k0=1)
        pl.plot_2D(np.angle(beam[:, :, -1]), map='Greys')
        pl.plot_2D(np.abs(beam[:, :, -1]))
        sing.Jz_calc_no_conj(beam[:, :, -1])
        # dots = sing.get_singularities(np.angle(beam))
        # pl.plot_2D(np.abs(beam[:, :, zRes//2]))
        # pl.plot_scatter_3D(dots[:, 0], dots[:, 1], dots[:, 2])
        # fig = pl.plot_3D_dots_go(dots)
        # pl.box_set_go(fig, mesh=None, autoDots=dots, perBox=0.05)
        # fig.show()
