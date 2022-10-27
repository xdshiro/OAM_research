import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np
def linear_momentum(EArray, xArray=None, yArray=None):
    EArray = np.array(EArray)
    Er, Ei = np.real(EArray), np.imag(EArray)
    if xArray is None or yArray is None:
        shape = np.shape(EArray)
        xArray = np.arange(shape[0])
        yArray = np.arange(shape[1])
    x0 = (xArray[-1] + xArray[0]) / 2
    y0 = (yArray[-1] + yArray[0]) / 2
    x = np.array(xArray) - x0
    y = np.array(yArray) - y0
    dx = xArray[1] - xArray[0]
    dy = yArray[1] - yArray[0]
    Px = np.zeros(np.shape(Er))
    Py = np.zeros(np.shape(Er))
    print(np.shape(Px))
    for i in range(1, len(xArray) - 1, 1):
        for j in range(1, len(yArray) - 1, 1):
            dErx = (Er[i + 1, j] - Er[i - 1, j]) / (2 * dx)
            dEry = (Er[i, j + 1] - Er[i, j - 1]) / (2 * dy)
            dEix = (Ei[i + 1, j] - Ei[i - 1, j]) / (2 * dx)
            dEiy = (Ei[i, j + 1] - Ei[i, j - 1]) / (2 * dy)
            Px[i, j] = Er[i, j] * dEix - Ei[i, j] * dErx
            Py[i, j] = Er[i, j] * dEiy - Ei[i, j] * dEry
    print(f'total Px={np.sum(Px)}, total Py={Py}')

    return Px, Py
if __name__ == '__main__':
    # w = 1.7
    # a00 = 1 - 2 * w ** 2 + 2 * w ** 4
    # a01 = 2 * w ** 2 - 4 * w ** 4
    # a02 = 2 * w ** 4
    # a20 = 4 * np.sqrt(2) * w ** 2
    # aCoeff = [a00, a01, a02, a20]
    # aSumSqr = 0.1 * np.sqrt(sum([a ** 2 for a in aCoeff]))
    # aCoeff /= aSumSqr
    # print(aCoeff)
    # print(f'a00 -> a01 -> a02 ->... -> a0n -> an0:')
    # exit()
    plot_trefoil = False
    if plot_trefoil:
        # xyzMesh = fg.create_mesh_XYZ(2.1, 2.1, 0.7, 50, 50, 50, zMin=None)
        xyzMesh = fg.create_mesh_XYZ(2.1, 2.1, 0.8, 90, 90, 90, zMin=None)
        beam = bp.LG_combination(*xyzMesh,
                                 coefficients=[1.71, -5.66, 6.38, -2.30, -4.36],
                                 # coefficients=[0,0,0,0, -4.36],
                                 modes=[(0, 0), (0, 1), (0, 2), (0, 3), (3, 0)],
                                 width=[1, 1, 1, 1, 1])

        cross = beam[:, :, 25]
        # pl.plot_2D(np.abs(cross))
        sing.Jz_calc_no_conj(cross)
        pl.plot_2D(np.abs(cross))
        pl.plot_2D(np.angle(cross))
        # Px, Py = linear_momentum(cross)
        # pl.plot_2D(Px)
        # pl.plot_2D(Py)
        # exit()
        dots = sing.get_singularities(np.angle(beam), axesAll=True)
        # pl.plot_2D(np.abs(beam[:, :, zRes//2]))
        # pl.plot_scatter_3D(dots[:, 0], dots[:, 1], dots[:, 2])
        fig = pl.plot_3D_dots_go(dots, marker={'size': 20, 'color': 'black',
                                               'line': dict(width=180, color='white')})
        pl.box_set_go(fig, mesh=None, autoDots=dots, perBox=0.05)
        fig.show()
    plot_hopf = False
    if plot_hopf:
        # xyzMesh = fg.create_mesh_XYZ(2.1, 2.1, 0.7, 50, 50, 50, zMin=None)
        xyzMesh = fg.create_mesh_XYZ(2.1, 2.1, 0.8, 90, 90, 90, zMin=None)
        beam = bp.LG_combination(*xyzMesh,
                                 coefficients=[3.12949298, -7.25104274,  4.38399864,  4.29059539],
                                 # coefficients=[0,0,0,0, -4.36],
                                 modes=[(0, 0), (0, 1), (0, 2), (2, 0)],
                                 width=[1, 1, 1, 1])

        # cross = beam[:, :, 25]
        # pl.plot_2D(np.abs(cross))
        # sing.Jz_calc_no_conj(cross)
        # pl.plot_2D(np.abs(cross))
        # pl.plot_2D(np.angle(cross))
        # Px, Py = linear_momentum(cross)
        # pl.plot_2D(Px)
        # pl.plot_2D(Py)
        # exit()
        dots = sing.get_singularities(np.angle(beam), axesAll=True)
        # pl.plot_2D(np.abs(beam[:, :, zRes//2]))
        # pl.plot_scatter_3D(dots[:, 0], dots[:, 1], dots[:, 2])
        fig = pl.plot_3D_dots_go(dots, marker={'size': 20, 'color': 'black',
                                               'line': dict(width=180, color='white')})
        pl.box_set_go(fig, mesh=None, autoDots=dots, perBox=0.05)
        fig.show()
    plot_unknot = True
    if plot_unknot:
        w = 2.5
        a1 = -1 + w ** 2
        a2 = -w ** 2
        a3 = -2 * w
        # xyzMesh = fg.create_mesh_XYZ(2.1, 2.1, 0.7, 50, 50, 50, zMin=None)
        xyzMesh = fg.create_mesh_XYZ(2.1, 2.1, 0.8, 90, 90, 90, zMin=None)
        beam = bp.LG_combination(*xyzMesh,
                                 coefficients=[a1, a2, a3],
                                 # coefficients=[0,0,0,0, -4.36],
                                 modes=[(0, 0), (0, 1), (1, 0)],
                                 width=[1, 1, 1])

        cross = beam[:, :, 25]
        # pl.plot_2D(np.abs(cross))
        sing.Jz_calc_no_conj(cross)
        pl.plot_2D(np.abs(cross))
        pl.plot_2D(np.angle(cross))
        # Px, Py = linear_momentum(cross)
        # pl.plot_2D(Px)
        # pl.plot_2D(Py)
        # exit()
        dots = sing.get_singularities(np.angle(beam), axesAll=True)
        # pl.plot_2D(np.abs(beam[:, :, zRes//2]))
        # pl.plot_scatter_3D(dots[:, 0], dots[:, 1], dots[:, 2])
        fig = pl.plot_3D_dots_go(dots, marker={'size': 20, 'color': 'black',
                                               'line': dict(width=180, color='white')})
        pl.box_set_go(fig, mesh=None, autoDots=dots, perBox=0.05)
        fig.show()