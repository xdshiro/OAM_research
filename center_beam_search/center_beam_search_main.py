import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np
from matplotlib.pyplot import polar


def LG_spectre_coeff(field, l, p, xM=(-1, 1), yM=(-1, 1), width=1., k0=1., mesh=None):
    if mesh is None:
        shape = np.shape(field)
        mesh = fg.create_mesh_XY(xMinMax=xM, yMinMax=yM, xRes=shape[0], yRes=shape[1])
        dS = ((xM[1] - xM[0]) / (shape[0] - 1)) * ((yM[1] - yM[0]) / (shape[1] - 1))
    else:
        xArray, yArray = fg.arrays_from_mesh(mesh)
        dS = (xArray[1] - xArray[0]) * (yArray[1] - yArray[0])
    # shape = np.shape(field)
    # xyMesh = fg.create_mesh_XY_old(xMax=xM[1], yMax=yM[1], xRes=shape[0], yRes=shape[1], xMin=xM[0], yMin=yM[0])
    LGlp = bp.LG_simple(*mesh, l=l, p=p, width=width, k0=k0)

    return np.sum(field * np.conj(LGlp)) * dS


def displacement_lateral(field_func, mesh, r_0, eta, **kwargs):
    xArray, yArray = fg.arrays_from_mesh(mesh)
    xArray, yArray = xArray - r_0 * np.cos(eta), yArray - r_0 * np.sin(eta)
    mesh_new = np.meshgrid(xArray, yArray, indexing='ij')
    return field_func(*mesh_new, **kwargs)


def displacement_deflection(field_func, mesh, eta, gamma, k=1, **kwargs):
    field_change = np.exp(1j * k * fg.rho(*mesh) * np.sin(gamma) * np.cos(fg.phi(*mesh) - eta))
    return field_func(*mesh, **kwargs) * field_change


def LG_transition_matrix_real_space(operator, l, p, l0, p0, xM=(-1, 1), yM=(-1, 1), shape=(100, 100),
                                    width=1., k0=1., mesh=None, **kwargs):
    if mesh is None:
        mesh = fg.create_mesh_XY(xM, yM, xRes=shape[0], yRes=shape[1])
        dS = ((xM[1] - xM[0]) / (shape[0] - 1)) * ((yM[1] - yM[0]) / (shape[1] - 1))
    else:
        xArray, yArray = fg.arrays_from_mesh(mesh)
        dS = (xArray[1] - xArray[0]) * (yArray[1] - yArray[0])
    operatorOnLG = operator(bp.LG_simple,
                            mesh, z=0, l=l0, p=p0, width=width, k0=k0, x0=0, y0=0, z0=0, **kwargs)
    LGlp = bp.LG_simple(*mesh, l=l, p=p, width=width, k0=k0)

    return np.sum(operatorOnLG * np.conj(LGlp)) * dS


# def LG_transition_matrix_real_space_fast(operator, l, p, l0, p0, xM=(-1, 1), yM=(-1, 1), shape=(100, 100),
#                                          mesh=None, **kwargs):
#     from scipy.special import jv
#     # if mesh is None:
#     #     mesh = fg.create_mesh_XY(xM, yM, xRes=shape[0], yRes=shape[1])
#     #     dS = ((xM[1] - xM[0]) / (shape[0] - 1)) * ((yM[1] - yM[0]) / (shape[1] - 1))
#     # else:
#     #     xArray, yArray = fg.arrays_from_mesh(mesh)
#     #     dS = (xArray[1] - xArray[0]) * (yArray[1] - yArray[0])
#     # operatorOnLG = operator(bp.LG_simple,
#     #                         mesh, z=0, l=l0, p=p0, width=1, k0=1, x0=0, y0=0, z0=0, **kwargs)
#     # LGlp = bp.LG_simple(*mesh, l=l, p=p, width=1, k0=1)
#
#     jv(2 * fg.rho(*mesh) * 1 / 1 ** 2)
#     return np.sum(operatorOnLG * np.conj(LGlp)) * dS

def variance_V_helper(Pl, lArray):
    sum1, sum2 = 0, 0
    for p, l in zip(Pl, lArray):
        sum1 += p * l ** 2
        sum2 += p * l
    return sum1 - sum2 ** 2


def variance_single_transition(field, mesh, r, eta, displacement_function=displacement_lateral,
                               p=(0, 5), l=(0, 4), p0=(0, 5), l0=(-3, 3),
                               width=1., k0=1.):
    p1, p2 = p
    l1, l2 = l
    p01, p02 = p0
    l01, l02 = l0
    Pl = np.zeros(l2 - l1 + 1)
    for ind_l, l in enumerate(np.arange(l1, l2 + 1)):
        sum = 0
        for p in np.arange(p1, p2 + 1):
            sum_inner = 0
            for l0 in np.arange(l01, l02 + 1):
                for p0 in np.arange(p01, p02 + 1):
                    if displacement_function is displacement_lateral:
                        element_ = LG_transition_matrix_real_space(displacement_function, l=l, p=p,
                                                                   l0=l0, p0=p0,
                                                                   mesh=mesh, r_0=r, eta=eta,
                                                                   width=width, k0=k0)
                    else:
                        element_ = LG_transition_matrix_real_space(displacement_function, l=l, p=p,
                                                                   l0=l0, p0=p0,
                                                                   mesh=mesh, eta=r, gamma=eta,
                                                                   width=width, k0=k0)
                    value = LG_spectre_coeff(field, l=l0, p=p0, mesh=mesh, width=width, k0=k0)
                    # print(element_, value, (np.abs(element_ * value) ** 2))
                    sum_inner += (element_ * value)
            sum += np.abs(sum_inner) ** 2
        Pl[ind_l] = sum
        # print(Pl[ind_l])
    V = variance_V_helper(Pl, np.arange(l1, l2 + 1))
    return V


# def variance_single_transition_tilt(field, mesh, r, eta, displacement_function=displacement_lateral,
#                                      p=(0, 5), l=(0, 4), p0=(0, 5), l0=(-3, 3),
#                                      width=1., k0=1.):
#
#     p1, p2 = p
#     l1, l2 = l
#     p01, p02 = p0
#     l01, l02 = l0
#     Pl = np.zeros(l2 - l1 + 1)
#     for ind_l, l in enumerate(np.arange(l1, l2 + 1)):
#         sum = 0
#         for p in np.arange(p1, p2 + 1):
#             sum_inner = 0
#             for l0 in np.arange(l01, l02 + 1):
#                 for p0 in np.arange(p01, p02 + 1):
#                     element_ = LG_transition_matrix_real_space(displacement_function, l=l, p=p,
#                                                                l0=l0, p0=p0,
#                                                                mesh=mesh, r_0=r, eta=eta,
#                                                                width=width, k0=k0)
#                     value = LG_spectre_coeff(field, l=l0, p=p0, mesh=mesh, width=width, k0=k0)
#                     # print(element_, value, (np.abs(element_ * value) ** 2))
#                     sum_inner += (element_ * value)
#             sum += np.abs(sum_inner) ** 2
#         Pl[ind_l] = sum
#         # print(Pl[ind_l])
#     V = variance_V_helper(Pl, np.arange(l1, l2 + 1))
#     return V


def variance_map_shift(beam, mesh, displacement_function=displacement_lateral,
                       resolution_V=(4, 4), xBound=(-1, 1), yBound=(-1, 1), width=1., k0=1.,
                       p=(0, 5), l=(0, 4), p0=(0, 5), l0=(-3, 3)):
    V = np.zeros(resolution_V)

    xArray = np.linspace(*xBound, resolution_V[0])
    yArray = np.linspace(*yBound, resolution_V[1])
    for i, x in enumerate(xArray):
        print('Main Coordinate x: ', i)
        for j, y in enumerate(yArray):
            print('y: ', j)
            r = fg.rho(x, y)
            eta = np.angle(x + 1j * y)
            V[i, j] = variance_single_transition(beam, mesh, r, eta,
                                                 displacement_function=displacement_function,
                                                 width=width, k0=k0,
                                                 p=p, l=l, p0=p0, l0=l0)
    return V


def variance_map_tilt(beam, mesh, displacement_function=displacement_deflection,
                      resolution_V=(4, 4), etaBound=(-1, 1), gammaBound=(-1, 1), width=1., k0=1.,
                      p=(0, 5), l=(0, 4), p0=(0, 5), l0=(-3, 3)):
    V = np.zeros(resolution_V)

    etaArray = np.linspace(*etaBound, resolution_V[0])
    gammaArray = np.linspace(*gammaBound, resolution_V[1])
    for i, eta in enumerate(etaArray):
        print('Main Coordinate eta: ', i)
        for j, gamma in enumerate(gammaArray):
            print('gamma: ', j)
            V[i, j] = variance_single_transition(beam, mesh, eta, gamma,
                                                 displacement_function=displacement_function,
                                                 width=width, k0=k0,
                                                 p=p, l=l, p0=p0, l0=l0)
    return V


def LG_spectrum(beam, l=(-3, 3), p=(0, 5), xM=(-1, 1), yM=(-1, 1), width=1., k0=1., mesh=None, plot=True):
    l1, l2 = l
    p1, p2 = p
    spectrum = np.zeros((l2 - l1 + 1, p2 - p1 + 1), dtype=complex)
    # spectrumReal = []
    # modes = []
    for l in np.arange(l1, l2 + 1):
        for p in np.arange(p1, p2 + 1):
            value = LG_spectre_coeff(beam, l=l, p=p, xM=xM, yM=yM, width=width, k0=k0, mesh=mesh)
            # print(l, p, ': ', value, np.abs(value))
            spectrum[l - l1, p] = value
            # if np.abs(value) > 0.5:
            # spectrumReal.append(value)
            # modes.append((l, p))
    # print(modes)
    if plot:
        pl.plot_2D(np.abs(spectrum), x=np.arange(l1 - 0.5, l2 + 1 + 0.5), y=np.arange(p1 - 0.5, p2 + 1 + 0.5),
                   interpolation='none', grid=True, xname='l', yname='p', show=False)
        plt.yticks(np.arange(p1, p2 + 1))
        plt.xticks(np.arange(l1, l2 + 1))
        plt.show()
    return spectrum


def center_beam_finding(beam, mesh, stepXY=None, displacement_function=displacement_lateral,
                        p=(0, 5), l=(0, 4), p0=(0, 5), l0=(-3, 3),
                        width=1, k0=1, x=None, y=None):
    if stepXY is None:
        xArray, yArray = fg.arrays_from_mesh(mesh)
        # print(xArray)
        stepXY = xArray[1] - xArray[0], yArray[1] - yArray[0]
    if x is None:
        xArray, yArray = fg.arrays_from_mesh(mesh)
        x = xArray[len(xArray) // 2]
    if y is None:
        xArray, yArray = fg.arrays_from_mesh(mesh)
        y = yArray[len(yArray) // 2]

    def search(xFlag):
        nonlocal x, y, var0
        varIt = var0
        signX = 1
        correctWay = False
        while True:
            if xFlag:
                x = x + signX * stepXY[0]
            else:
                y = y + signX * stepXY[1]
            var = variance_single_transition(beam, mesh=xyMesh, displacement_function=displacement_function,
                                             r=fg.rho(x, y), eta=np.angle(x + 1j * y),
                                             p=p, l=l, p0=p0, l0=l0, width=width, k0=k0)
            print(x, y, var)
            if var < varIt:
                varIt = var
                correctWay = True
            else:
                if xFlag:
                    x = x - signX * stepXY[0]
                else:
                    y = y - signX * stepXY[1]
                if correctWay:
                    break
                else:
                    correctWay = True
                    signX *= -1

        return varIt

    var0 = variance_single_transition(beam, mesh=xyMesh, r=fg.rho(x, y), eta=np.angle(x + 1j * y),
                                      p=p, l=l, p0=p0, l0=l0, width=width, k0=k0)
    print(f'var0={var0}')

    while True:
        search(xFlag=True)
        varXY = search(xFlag=False)
        if varXY == var0:
            return -x, -y
        else:
            var0 = varXY


def tilt_beam_finding(beam, mesh, stepEG=None, displacement_function=displacement_deflection,
                        p=(0, 5), l=(0, 4), p0=(0, 5), l0=(-3, 3),
                        width=1, k0=1, eta=0, gamma=0):
    if stepEG is None:
        stepEG = 1, 1

    def search(etaFlag):
        nonlocal eta, gamma, var0
        varIt = var0
        signX = 1
        correctWay = False
        while True:
            if etaFlag:
                eta = eta + signX * stepEG[0]
            else:
                gamma = gamma + signX * stepEG[1]
            var = variance_single_transition(beam, mesh, eta, gamma,
                                             displacement_function=displacement_function,
                                             p=p, l=l, p0=p0, l0=l0, width=width, k0=k0)

            print(eta, gamma, var)
            if var < varIt:
                varIt = var
                correctWay = True
            else:
                if etaFlag:
                    eta = eta - signX * stepEG[0]
                else:
                    gamma = gamma - signX * stepEG[1]
                if correctWay:
                    break
                else:
                    correctWay = True
                    signX *= -1

        return varIt

    var0 = variance_single_transition(beam, mesh, eta, gamma,
                                             displacement_function=displacement_function,
                                             p=p, l=l, p0=p0, l0=l0, width=width, k0=k0)
    print(f'var0={var0}')

    while True:
        search(etaFlag=True)
        varEG = search(etaFlag=False)
        if varEG == var0:
            return -eta, -gamma
        else:
            var0 = varEG


def find_width(beam, mesh, widthStep=0.1, l=(-3, 3), p=(0, 5), width=1., k0=1.):
    """
    this function finds the approximate beam waste (any position, any beam=sum(LG))
    :param mesh: mesh is required to have a scale of the beam => correct width scale
    :param widthStep: precision of the beam width search (but it's not real precision, method is not very accurate)
    :param l: we wont to cover all spectrum
    :param p: we wont to cover all spectrum
    :param width: starting beam width (you better have it as accurate as poosible, since the method is bad)
    :param k0: does nothing
    :return: ~beam width
    """
    minSpec = np.sum(np.abs(LG_spectrum(beam, l=l, p=p, mesh=mesh, width=width, k0=k0, plot=False)) ** (1 / 2))
    correctWay = False
    direction = +1
    while True:
        width += direction * widthStep
        spec = np.sum(np.abs(LG_spectrum(beam, l=l, p=p, mesh=mesh, width=width, k0=k0, plot=False)) ** (1 / 2))
        print(width, spec)
        if spec < minSpec:
            minSpec = spec
            correctWay = True
        else:
            width -= direction * widthStep
            if correctWay:
                break
            else:
                correctWay = True
                direction *= -1
    return width


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    xB, yB = [-4, 4], [-4, 4]
    # xyzMesh = fg.create_mesh_XYZ(2.1, 2.1, 0.7, 50, 50, 50, zMin=None)
    xyMesh = fg.create_mesh_XY_old(xB[1], yB[1], 50, 50, xMin=xB[0], yMin=yB[0])


    def beamF(*xyMesh, **kwargs):
        return bp.LG_combination(*xyMesh,
                                 coefficients=[1, 1],
                                 modes=[(0, 0), (2, 1)],
                                 width=[1, 1])


    xDis = 0.5
    yDis = 0.5
    eta = -50
    gamma = 20  # degrees
    beam_displaced = displacement_lateral(beamF, xyMesh, r_0=fg.rho(xDis, yDis),
                                          eta=np.angle(xDis + 1j * yDis))
    beam = beam_displaced
    beam_tilted = displacement_deflection(beamF, xyMesh,
                                          eta=eta * np.pi / 180, gamma=gamma * np.pi / 180)
    # beam = beam_tilted
    # pl.plot_2D(np.angle(beam), x=xB, y=yB)
    pl_dict = {'p': (0, 3), 'l': (-3, 4), 'p0': (0, 3), 'l0': (-3, 4)}

    width = 1.0
    k0 = 1
    # width_close = find_width(beam, mesh=xyMesh, widthStep=0.03, l=(-7, 7), p=(0, 7), width=0.5, k0=1.)
    # print(width_close)
    # exit()
    # k0 = 1
    # spec = LG_spectrum(beam, l=(-3, 4), p=(0, 3), mesh=xyMesh, plot=True, width=width, k0=k0)
    # exit()
    # V = variance_map_tilt(beam=beam, mesh=xyMesh,
    #                       resolution_V=(20, 20), etaBound=(30 * np.pi / 180, 60 * np.pi / 180),
    #                       gammaBound=(-30 * np.pi / 180, 0 * np.pi / 180),
    #                       **pl_dict, width=width)
    # print(V)
    # pl.plot_2D(V, x=[30, 60], y=[-30,0])
    # exit()
    # exit()
    # print(np.sum(np.abs(spec)))
    # exit()
    # V = variance_map_shift(beam=beam, mesh=xyMesh,
    #                  resolution_V=(5, 5), xBound=(-1, 1), yBound=(-1, 1),
    #                  **pl_dict, width=width)
    # print(V)
    # pl.plot_2D(V, x=[-1, 1], y=[-1, 1])
    # exit()
    #
    # x = -0.5
    # y = -0.5
    # print(variance_single_transition(beam, mesh=xyMesh,
    #                                  r=fg.rho(x, y), eta=np.angle(x + 1j * y),
    #                                  **pl_dict,
    #                                  width=1, k0=1))
    # x = -0.4
    # y = -0.4
    # print(variance_single_transition(beam, mesh=xyMesh,
    #                                  r=fg.rho(x, y), eta=np.angle(x + 1j * y),
    #                                  **pl_dict,
    #                                  width=1, k0=1))
    #
    # exit()
    x, y = center_beam_finding(beam, xyMesh, x=0, y=0, **pl_dict, stepXY=[0.1, 0.1])
    print(x, y)
    print(eta * np.pi / 180, gamma * np.pi / 180)
    print(tilt_beam_finding(beam, xyMesh, eta=-30 * np.pi / 180, gamma=-15 * np.pi / 180, **pl_dict, stepEG=[0.05, 0.05]))
    exit()

    check_spectrum = 1
    if check_spectrum:
        l1, l2 = -3, 3
        p1, p2 = 0, 3
        spectrum = np.zeros((l2 - l1 + 1, p2 - p1 + 1))
        spectrumReal = []
        modes = []
        for l in np.arange(l1, l2 + 1):
            for p in np.arange(p1, p2 + 1):
                value = LG_spectre_coeff(beam, l=l, p=p, xM=xB, yM=yB, width=1, k0=1)
                # print(l, p, ': ', value, np.abs(value))
                spectrum[l - l1, p] = np.abs(value)
                # if np.abs(value) > 0.5:
                spectrumReal.append(value)
                modes.append((l, p))
        # print(modes)

        pl.plot_2D(spectrum, y=np.arange(l1 - 0.5, l2 + 1 + 0.5), x=np.arange(p1 - 0.5, p2 + 1 + 0.5),
                   interpolation='none', grid=True, xname='p', yname='l', show=False)
        plt.yticks(np.arange(p1, p2 + 1))
        plt.xticks(np.arange(l1, l2 + 1))
        plt.show()
        exit()
        beam = bp.LG_combination(*xyMesh, coefficients=spectrumReal, modes=modes)
        pl.plot_2D(np.abs(beam), axis_equal=True)
        pl.plot_2D(np.angle(beam), axis_equal=True, vmin=-np.pi, vmax=np.pi)
        sing.Jz_calc_no_conj(beam)

    # beam_displaced = displacement_lateral(gauss, xyMesh, r_0=1, eta=0.3)
    # pl.plot_2D(np.abs(beam_displaced), x=[-3, 3], y=[-3, 3], show=True, axis_equal=True)
    # beam_deflected = displacement_deflection(beam, xyMesh, gamma=0.5, eta=0.4)
    # pl.plot_2D(np.angle(beam_deflected), x=xB, y=yB, show=True, axis_equal=True)
    exit()

    plot_trefoil = True
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
