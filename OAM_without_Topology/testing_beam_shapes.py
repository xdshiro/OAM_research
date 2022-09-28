import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np


def OAM_no_topology_test(x, y, z=0, l=1, p=0, width=1, k0=1, x0=0, y0=0, z0=0):
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

    x = x - x0
    y = y - y0
    z = z - z0

    def f2(phi):
        return np.abs(phi)

    def f1(phi, A0=2, A1=1, A2=1):
        phi[phi >= 0] = A0 + A1 * np.sin(phi[phi >= 0])
        phi[phi < 0] = A0 + A2 * np.sin(phi[phi < 0])
        return phi

    E = f1(fg.phi(x, y)) * np.exp(1j * f2(fg.phi(x, y))) * np.exp(-fg.rho(x, y) ** 2 / (2 * width ** 2))
    # zR = k0 * width ** 2
    # E = (np.sqrt(np.math.factorial(p) / (np.pi * np.math.factorial(np.abs(l) + p)))
    #      * fg.rho(x, y) ** np.abs(l) * np.exp(1j * l * fg.phi(x, y))
    #      / (width ** (np.abs(l) + 1) * (1 + 1j * z / zR) ** (np.abs(l) + 1))
    #      * ((1 - 1j * z / zR) / (1 + 1j * z / zR)) ** p
    #      * np.exp(-fg.rho(x, y) ** 2 / (2 * width ** 2 * (1 + 1j * z / zR)))
    #      * laguerre_polynomial(fg.rho(x, y) ** 2 / (width ** 2 * (1 + z ** 2 / zR ** 2)), np.abs(l), p)
    #      )
    return E


def LG_spectr_coeff(field, l, p, xM=(-1, 1), yM=(-1, 1), width=1, k0=1):
    shape = np.shape(field)
    xyMesh = fg.create_mesh_XY_old(xMax=xM[1], yMax=yM[1], xRes=shape[0], yRes=shape[1], xMin=xM[0], yMin=yM[0])
    LGlp = bp.LG_simple(*xyMesh, l=l, p=p, width=width, k0=k0)
    dS = ((xM[1] - xM[0]) / (shape[0] - 1)) * ((yM[1] - yM[0]) / (shape[1] - 1))
    return np.sum(field * np.conj(LGlp)) * dS


if __name__ == '__main__':
    test_spectrum = True
    if test_spectrum:
        xyMesh = fg.create_mesh_XY_old(12.5, 12.5, 350, 350)
        # beam = bp.LG_simple(*xyMesh, l=1, p=1, width=1, k0=1) * np.exp(1j * 0.4)
        beam = OAM_no_topology_test(*xyMesh)
        beam = fg.propagator_split_step_3D(beam, dz=1, xArray=None, yArray=None,
                                           zSteps=300, n0=1, k0=1)[:, :, -1]
        pl.plot_2D(np.abs(beam))
        pl.plot_2D(np.angle(beam))
        sing.Jz_calc_no_conj(beam)
        exit()
        l1, l2 = -5, 5
        p1, p2 = 0, 30
        z = np.zeros((l2 - l1 + 1, p2 - p1 + 1))
        zReal = []
        modes = []
        for l in np.arange(l1, l2 + 1):
            for p in np.arange(p1, p2 + 1):
                value = LG_spectr_coeff(beam, l=l, p=p, xM=(-7.5, 7.5), yM=(-7.5, 7.5), width=1, k0=1)
                # print(l, p, ': ', value, np.abs(value))
                z[l- l1,p]=np.abs(value)
                # if np.abs(value) > 0.5:
                zReal.append((value))
                modes.append((l, p))
        # print(modes)
        pl.plot_2D(z, y=np.arange(l1, l2 + 1), x= np.arange(p1, p2 + 1), interpolation='none')
        beam = bp.LG_combination(*xyMesh,coefficients=zReal, modes=modes)
        beam = fg.propagator_split_step_3D(beam, dz=0.5, xArray=None, yArray=None,
                                           zSteps=300, n0=1, k0=1)[:, :, -1]
        # pl.plot_2D(np.abs(beam)[125-30: 125+30, 125-30: 125+30])
        # pl.plot_2D(np.angle(beam)[125-30: 125+30, 125-30: 125+30], vmin=-np.pi, vmax=np.pi)
        # sing.Jz_calc_no_conj(beam[125-30: 125+30, 125-30: 125+30])
        pl.plot_2D(np.abs(beam))
        pl.plot_2D(np.angle(beam), vmin=-np.pi, vmax=np.pi)
        sing.Jz_calc_no_conj(beam)
        # pl.plot_2D(np.abs(beam))
        # beam = fg.propagator_split_step_3D(beam, dz=5.25,
        #                                    zSteps=70, n0=1, k0=1)
        # beam = beam[:, :, -1]
        width = 1
        # width = np.sqrt(width ** 2 * (1 + (0.25 * 70) ** 2 / 1 ** 2))
        # pl.plot_2D(np.abs(beam))
        # print(width)
        # print(LG_spectr_coeff(beam, l=1, p=1, xM=(-7.5, 7.5), yM=(-7.5, 7.5), width=width, k0=1))


    copy_to_another_file = False
    if copy_to_another_file:
        xyMesh = fg.create_mesh_XY_old(7.5, 7.5, 350, 350)
        beam = OAM_no_topology_test(*xyMesh)
        # pl.plot_2D(np.abs(beam))
        # pl.plot_2D(np.angle(beam))
        beam = fg.propagator_split_step_3D(beam, dz=0.25, xArray=None, yArray=None, zSteps=70, n0=1, k0=1)
        # pl.plot_2D(np.angle(beam[:, :, -1]), map='jet')
        # pl.plot_2D(np.abs(beam[:, :, -1]))
        cr = 60
        new_beam = fg.cut_filter(beam[350 // 2 - cr:350 // 2 + cr, 350 // 2 - cr:350 // 2 + cr, -1], radiusPix=60,
                                 circle=False, phaseOnly=False)
        # pl.plot_2D(np.angle(new_beam), map='jet')
        # pl.plot_2D(np.abs(new_beam))
        sing.Jz_calc_no_conj(beam[:, :, -1])
    plot_trefoil = False
    if plot_trefoil:
        xyzMesh = fg.create_mesh_XYZ(2.1, 2.1, 0.7, 50, 50, 50, zMin=None)
        beam = bp.LG_combination(*xyzMesh,
                                 coefficients=[1.71, -5.66, 6.38, -2.30, -4.36],
                                 modes=[(0, 0), (0, 1), (0, 2), (0, 3), (3, 0)],
                                 width=[1, 1, 1, 1, 1])
        dots = sing.get_singularities(np.angle(beam))
        # pl.plot_2D(np.abs(beam[:, :, zRes//2]))
        # pl.plot_scatter_3D(dots[:, 0], dots[:, 1], dots[:, 2])
        fig = pl.plot_3D_dots_go(dots)
        pl.box_set_go(fig, mesh=None, autoDots=dots, perBox=0.05)
        fig.show()
