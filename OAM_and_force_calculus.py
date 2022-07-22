import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go


# fig = go.Figure(data=go.Cone(x=[1], y=[1], z=[1], u=[1], v=[1], w=[0]))
#
# fig.update_layout(scene_camera_eye=dict(x=-0.76, y=1.8, z=0.92))
def plot_3D_cones(data, mesh, colorscale='Blues', indexing='ij',
                  sizemode="absolute", sizeref=1., filterValue=None, power=1, random=1):
    if indexing != 'ij':
        print(f'no "xy" version')
        exit()

    def _transposition_ij_xy(md):
        shapeOld = np.shape(md)
        shape = (shapeOld[0], shapeOld[2], shapeOld[1], shapeOld[3])
        ans = np.zeros((shape))
        for i in range(shapeOld[0]):
            for j in range(shapeOld[3]):
                ans[i][:, :, j] = md[i][:, :, j].transpose()
        return ans

    def _narray_flattening(mesh):
        ans = []
        for m in mesh:
            ans.append(m.flatten())
        return ans

    def _filter_data_for_cones(data, mesh, valueMin):
        d = _narray_flattening(_transposition_ij_xy(data))
        m = _narray_flattening(_transposition_ij_xy(mesh))

        dataAbs = np.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2)
        fullData = np.array(list(zip(dataAbs / dataAbs.max(), d[0], d[1], d[2], m[0], m[1], m[2])))
        if random < 1 and random > 0:
            elements = list(range(len(fullData)))
            elementsRandom = np.random.choice(elements, size=int(len(fullData) * (1 - random)), replace=False)
            fullData = np.delete(fullData, elementsRandom, axis=0)

        elementsDelete = []
        for i, elem in enumerate(fullData):
            if elem[0] < valueMin:
                elementsDelete.append(i)

        filteredData = np.delete(fullData, elementsDelete, axis=0)
        filteredData[:, 1:4] /= dataAbs.max()
        return ((filteredData[:, 1], filteredData[:, 2], filteredData[:, 3]),
                (filteredData[:, 4], filteredData[:, 5], filteredData[:, 6]))

    if isinstance(data, list):
        data = np.array(data)
    if filterValue is not None or (random < 1 and random > 0):
        d, m = _filter_data_for_cones(data, mesh, valueMin=filterValue)
    else:
        d = _narray_flattening(_transposition_ij_xy(data))
        m = _narray_flattening(_transposition_ij_xy(mesh))
    if power != 1:
        dataAbs = np.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2) ** power
        d = [dd * dataAbs for dd in d]
    fig = go.Figure(data=go.Cone(
        x=m[0], y=m[1], z=m[2],
        u=d[0], v=d[1], w=d[2],
        colorscale=colorscale,
        sizemode=sizemode,
        sizeref=sizeref,
        opacity=1,
        anchor="tip"))
    fig.update_layout(
        scene=dict(domain_x=[0, 1],
                   camera_eye=dict(x=-1.57, y=1.36, z=0.58)))

    return fig

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


def poynting_vector_avg(field):
    # 0.5 * np.real(E x np.conj(H))
    # E cross H / 2
    pass


def force_optical(field, mesh, alpha=0 + 0j, c=1, indexing='ij'):
    # DOI: 10.1103/PhysRevApplied.14.034069
    # DOI: 10.1103/PhysRevLett.102.113602
    xyz = fg.arrays_from_mesh(mesh, indexing=indexing)

    # xyz = mesh
    def gradient_force(field, alpha):
        """
        gradient force is proportional to the gradient of the potential energy of the particle under the influence
         of the electromagnetic field (gradI)
        """
        return np.gradient(np.real(alpha) / 4 * np.abs(field) ** 2, *xyz)

    def scattering_force(field, c):
        """
        scattering force and is proportional to the Poynting vector (for Gaussian beam the scattering force
        points in the direction of propagation of the beam near the origin)
        """
        # sigma * 1 / c * poynting_vector(field)
        # fake
        # np.imag(1 * np.conj(fieldOLD[:, yPar, :]) * np.gradient(fieldOLD[:, yPar, :]))
        grad = np.gradient(field, *xyz)
        ans = [np.imag(alpha) / 2  * np.imag(np.conj(field) * g) for g in grad]
        # force = [np.imag(sigma * field * g) for g in grad]
        return ans

    def spin_curl_force():
        """
        spin-curl force, is a result of polarization gradients and can be disregarded in the case of uniform
        linear polarization in which we are interested.
        """
        return 0

    # force = (gradient_force(field, alpha)[0] + scattering_force(field, sigma, c)[0],
    #          gradient_force(field, alpha)[1] + scattering_force(field, sigma, c)[1])
    forceGrad = gradient_force(field, alpha)
    forceScat = scattering_force(field, c)
    forceTotal = tuple((forceGrad[i] + forceScat[i]) for i in range(len(forceGrad)))
    return forceTotal


def plot_vector_field_2D(*vectorField, mesh=None, show=True, density=1, color=None,
                         lw=None, cmap=None, indexing='ij'):
    """
    function plots 2D vector field
    :param vectorField: Fx, Fy
    :param xy: Mesh in 'ij' indexing
    :param density: density of the lines
    """
    if indexing != 'ij':
        print(f'no "xy" version')
        exit()
    if mesh is None:
        shape = np.shape(vectorField[0])
        mesh = fg.create_mesh_XY(1, 1, shape[0], shape[1], xMin=0, yMin=0, indexing=indexing)

    absValue = np.sqrt(vectorField[0] ** 2 + vectorField[1] ** 2)
    print(f'max force = {absValue.max()}')
    if lw is None:
        lw = (absValue / absValue.max() * 4).transpose()
    if color is None:
        color = (absValue / absValue.max() * 4).transpose()
    # xy = fg.arrays_from_mesh(mesh, indexing=indexing)
    xy = mesh
    # print(np.shape(vectorField[0]), np.shape(lw))
    # strm = plt.streamplot(xy[1], xy[0], vectorField[0].transpose(), vectorField[1].transpose(),
    #                       density=density, color=color, linewidth=lw, cmap=cmap)
    strm = plt.streamplot(xy[0].transpose(), xy[1].transpose(),
                          vectorField[0].transpose(), vectorField[1].transpose(),  # .transpose()
                          density=density, color=color, linewidth=lw, cmap=cmap)

    if cmap is not None:
        plt.colorbar(strm.lines)
    if show:
        plt.show()


# xMax = 5
# xMin = -3
# pl = ((xMax+xMin) / 2 +  (np.random.rand(200) - 0.5) * (xMax - xMin))
# plt.scatter(np.linspace(xMin, xMax, 200), pl)
# plt.show()
# exit()

if __name__ == '__main__':
    # a = [[1, 2], [1,3]]
    # b = [[1, 2], [1,3]]
    # print(np.cross(a))
    # exit()
    xMinMax = 3.5
    yMinMax = 3.5
    zMinMax = 0.5
    zRes = 140  # 250
    xRes = 140
    yRes = 140

    indexing = 'ij'
    xyzMesh = fg.create_mesh_XYZ(xMinMax, yMinMax, zMinMax, xRes, yRes, zRes, indexing=indexing)
    xyMesh = fg.create_mesh_XY(xMinMax, yMinMax, xRes, yRes, indexing=indexing)
    # print(xyzMesh)
    # modes = [(0, 1), (1, 0)]
    # coeff = [0, 1]
    modes = ((0, 0), (0, 1), (0, 2), (0, 3), (3, 0))
    coeff = [1.715, -5.662, 6.381, -2.305, -4.356]
    # coeff = [0, 0, 0, 0, 1]
    # beam = bp.LG_combination(*xyMesh, coefficients=coeff, modes=modes)
    beam = bp.LG_combination(*xyzMesh, coefficients=coeff, modes=modes)
    force = force_optical(beam, mesh=xyzMesh, alpha=-1 + 0.0001j, indexing=indexing)
    #
    #
    # def _transposition_ij_xy(md):
    #     shapeOld = np.shape(md)
    #     shape = (shapeOld[0], shapeOld[2], shapeOld[1], shapeOld[3])
    #     ans = np.zeros((shape))
    #     for i in range(shapeOld[0]):
    #         for j in range(shapeOld[3]):
    #             ans[i][:, :, j] = md[i][:, :, j].transpose()
    #     return ans
    #

    # mesh = _transposition_ij_xy(xyzMesh)
    # force = _transposition_ij_xy(force)
    # xyzMesh = fg.create_mesh_XYZ(xMinMax, yMinMax, zMinMax, xRes, yRes, zRes, indexing='xy')
    # pl.plot_2D(abs(beam[:, :, zRes // 2]), show=True)
    fig = plot_3D_cones(force, xyzMesh, sizeref=14, filterValue=0.1, colorscale='jet', power=2, random=0.002)
    dots = sing.get_singularities(np.angle(beam), bigSingularity=False, axesAll=True)
    dots = fg.dots3D_rescale(dots, xyzMesh)
    fig.add_trace(pl.plot_3D_dots_go(dots))
    fig.show()
    # plot_vector_field_2D(force[0][:, :, zRes // 2], force[1][:, :, zRes // 2], lw=2,
    #                      indexing=indexing, density=2, show=True, color=None, cmap='jet')
    # plot_vector_field_2D(force[0][:, :, zRes // 2], force[1][:, :, zRes // 2],
    #                      mesh=(mesh[0][:, :, zRes // 2], mesh[1][:, :, zRes // 2]),
    #                      indexing=indexing, density=3, show=True, color='k', cmap='BuPu')

    exit()
    # xyz = fg.arrays_from_mesh(xyzMesh)
    # beam = bp.LG_combination(*xyz, coefficients=coeff, modes=modes)
    # print(xy)
    # plt.scatter(xy[0],abs(beam[:, 19, 19]))
    # plt.show()
    # exit()
    # pl.plot_2D(np.angle(beam))
    # print(fg.arrays_from_mesh(xyzMesh))
    # exit()
    # pl.plot_scatter_3D(*xyz)
    # force = force_optical(beam, mesh=xyzMesh, alpha=1)
    # plot_vector_field_2D(*force, density=2, show=True,
    #                      color=np.sqrt(np.abs(force[0]) ** 2 + np.abs(force[1]) ** 2), cmap='BuPu')
    # print(grad)
    # beam = bp.LG_combination(*xyMesh, coefficients=coeff, modes=modes)
    # xy = fg.arrays_from_mesh(xyMesh)
    # beam = [[1,3,5],[2,6,9],[3,10,19]]
    # pl.plot_2D(beam)
    # pl.plot_2D(force[1])

    exit()
    plot_vector_field_2D(*force, mesh=xyzMesh, density=3, show=True, indexing=indexing,
                         color=np.sqrt(np.abs(force[0]) ** 2 + np.abs(force[1]) ** 2), cmap='BuPu')
    exit(0)
