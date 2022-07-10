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
def plot_3D_cones(data, mesh, colorscale='Blues', sizemode="absolute", sizeref=1., filterValue=None):
    def _narray_flattening(mesh):
        ans = []
        for m in mesh:
            ans.append(m.flatten())
        return ans

    def _filter_data_for_cones(data, mesh, valueMin):
        d = _narray_flattening(data)
        m = _narray_flattening(mesh)
        dataAbs = np.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2)
        fullData = np.array(list(zip(dataAbs / dataAbs.max(), d[0], d[1], d[2], m[0], m[1], m[2])))
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
    d, m = _filter_data_for_cones(data, mesh, valueMin=0.3)
    fig = go.Figure(data=go.Cone(
        x=m[0], y=m[1], z=m[2],
        u=d[0], v=d[1], w=d[2],
        colorscale=colorscale,
        sizemode=sizemode,
        sizeref=sizeref,
        anchor="tip"))
    fig.update_layout(
        scene=dict(domain_x=[0, 1],
                   camera_eye=dict(x=-1.57, y=1.36, z=0.58)))

    fig.show()
    # x = np.linspace(-1, 1, 15)
    # y = np.linspace(-1, 1, 15)
    # z = np.linspace(0.1, 0.3, 15)
    # xyzMesh = np.meshgrid(x,y,z, indexing='ij')
    # # noinspection PyTypeChecker
    # values = abs(xyzMesh[0]**2-xyzMesh[1]**2)/xyzMesh[2]
    # values = values/abs(values).max()
    # grad = np.gradient(values)
    #
    # exit()


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


def force_optical(field, alpha=0, sigma=0, c=1):
    # DOI: 10.1103/PhysRevApplied.14.034069
    # DOI: 10.1103/PhysRevLett.102.113602
    def gradient_force(field, alpha):
        """
        gradient force is proportional to the gradient of the potential energy of the particle under the influence
         of the electromagnetic field (gradI)
        """
        return np.gradient(np.real(alpha) / 4 * np.abs(field) ** 2)

    def scattering_force(field, sigma, c):
        """
        scattering force and is proportional to the Poynting vector (for Gaussian beam the scattering force
        points in the direction of propagation of the beam near the origin)
        """
        # sigma * 1 / c * poynting_vector(field)
        # fake
        # np.imag(1 * np.conj(fieldOLD[:, yPar, :]) * np.gradient(fieldOLD[:, yPar, :]))
        grad = np.gradient(np.conj(field))
        return np.imag(sigma * field * grad[0]), np.imag(sigma * field * grad[1])

    def spin_curl_force():
        """
        spin-curl force, is a result of polarization gradients and can be disregarded in the case of uniform
        linear polarization in which we are interested.
        """
        return 0

    # force = (gradient_force(field, alpha)[0] + scattering_force(field, sigma, c)[0],
    #          gradient_force(field, alpha)[1] + scattering_force(field, sigma, c)[1])
    force = gradient_force(field, alpha)
    return force


def plot_vector_field_2D(*vectorField, xy=None, show=True, density=1, color='k'):
    """
    function plots 2D vector field
    :param vectorField: Fx, Fy
    :param xy: Mesh in 'ij' indexing
    :param density: density of the lines
    """
    if xy is None:
        shape = np.shape(vectorField[0])
        xy = fg.create_mesh_XY(1, 1, shape[0], shape[1], xMin=0, yMin=0, indexing='ij')
    absValue = np.sqrt(vectorField[0] ** 2 + vectorField[1] ** 2)
    lw = absValue / absValue.max() * 6
    plt.streamplot(xy[1], xy[0], vectorField[0].transpose(), vectorField[1].transpose(),
                   density=density, color=color, linewidth=lw)
    if show:
        plt.show()


if __name__ == '__main__':
    # a = [[1, 2], [1,3]]
    # b = [[1, 2], [1,3]]
    # print(np.cross(a))
    # exit()
    xMinMax = 4
    yMinMax = 4
    zMinMax = 1
    zRes = 10
    xRes = yRes = 10
    xyzMesh = fg.create_mesh_XYZ(xMinMax, yMinMax, zMinMax, xRes, yRes, zRes, zMin=None)
    # xyMesh = fg.create_mesh_XY(xMinMax, yMinMax, xRes, yRes, indexing='ij')

    modes = [(0, 1), (1, 0)]
    coeff = [0, 1]
    beam = bp.LG_combination(*xyzMesh, coefficients=coeff, modes=modes)

    # pl.plot_2D(abs(beam))

    force = force_optical(beam, alpha=1, sigma=0)
    plot_3D_cones(force, xyzMesh, sizeref=1.5)
    # plot_vector_field_2D(*force)
    # print(grad)
    exit(0)
