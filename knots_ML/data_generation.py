"""
This script reads the field from mat file and produces all the necessary pre-processing procedures, and
creates a 3D array of singularity dots.

First main function is main_field_processing:
    1) reading the field from matlab file
    2) converting it into numpy array
    3) normalizing
    4) finding the beam waste
    5) rescaling the field, using the interpolation, for faster next steps
    6) finding the beam center
    7) rescaling field to the scale we want for 3D calculations
    8) removing the tilt and shift

Second main function is !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
"""
import my_functions.functions_general as fg
# import importlib
# importlib.reload(fg)
import my_functions.singularities as sing
import my_functions.beams_and_pulses as bp
import math
import numpy as np
import scipy.io as sio
import knots_ML.center_beam_search as cbs
import matplotlib.pyplot as plt


def read_field_2D_single(path, field=None):
    """
    Function reads .mat 2D array from matlab and convert it into numpy array

    If field is None, it will try to find the field name automatically

    :param path: full path to the file
    :param field: the name of the column with the field you want to read
    """
    field_read = sio.loadmat(path, appendmat=False)
    if field is None:
        for field_check in field_read:
            if len(np.shape(np.array(field_read[field_check]))) == 2:
                field = field_check
                break
    return np.array(field_read[field])


def normalization_field(field):
    """
    Normalization of the field for the beam center finding
    """
    field_norm = field / np.sqrt(np.sum(np.abs(field) ** 2))
    return field_norm


def plot_field(field):
    """
    Function plots intensity and phase of the field in 1 plot.
    Just a small convenient wrapper
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
    image1 = ax1.imshow(np.abs(field))
    ax1.set_title('|E|')
    plt.colorbar(image1, ax=ax1, shrink=0.4, pad=0.02, fraction=0.1)
    image2 = ax2.imshow(np.angle(field), cmap='jet')
    ax2.set_title('Phase(E)')
    plt.colorbar(image2, ax=ax2, shrink=0.4, pad=0.02, fraction=0.1)
    plt.tight_layout()
    plt.show()


def plot_field_3D_multi_planes(field3D, number=6, rows=3):
    """
    Function plots |E| and phase of the field in 1 plot.
    Just a small convenient wrapper
    """
    fig, axis = plt.subplots(math.ceil(number / rows), 3, figsize=(10, 3 * math.ceil(number / rows)))
    reso_z = np.shape(field3D)[2]
    for i, ax_r in enumerate(axis):
        for j, ax in enumerate(ax_r):
            image = ax.imshow(np.abs(field3D[:, :, int((reso_z - 1) / (number-1) * (len(ax_r) + j))]))
            ax.set_title(f'|E|, index z={int((reso_z - 1) / (number-1) * (i * len(ax_r) + j))}')
            plt.colorbar(image, ax=ax, shrink=0.4, pad=0.02, fraction=0.1)
    plt.tight_layout()
    plt.show()
    fig, axis = plt.subplots(math.ceil(number / rows), 3, figsize=(10, 3 * math.ceil(number / rows)))
    for i, ax_r in enumerate(axis):
        for j, ax in enumerate(ax_r):
            image = ax.imshow(np.angle(field3D[:, :, int((reso_z - 1) / (number - 1) * (len(ax_r) + j))]), cmap='jet')
            ax.set_title(f'phase(E), index z={int((reso_z - 1) / (number - 1) * (i * len(ax_r) + j))}')
            plt.colorbar(image, ax=ax, shrink=0.4, pad=0.02, fraction=0.1)

    plt.tight_layout()
    plt.show()


def find_beam_waist(field, mesh=None):
    """
    wrapper for the beam waste finder. More details in knots_ML.center_beam_search
    """
    shape = np.shape(field)
    if mesh is None:
        mesh = fg.create_mesh_XY(xRes=shape[0], yRes=shape[1])
    width = cbs.find_width(field, mesh=mesh, width=shape[1] // 8, widthStep=1, print_steps=False)
    return width


def field_interpolation(field, mesh=None, resolution=(100, 100),
                        xMinMax_frac=(1, 1), yMinMax_frac=(1, 1)):
    """
    Wrapper for the field interpolation fg.interpolation_complex
    :param resolution: new field resolution
    :param xMinMax_frac: new dimension for the field. (x_dim_old * frac)
    :param yMinMax_frac: new dimension for the field. (y_dim_old * frac)
    """
    shape = np.shape(field)
    if mesh is None:
        mesh = fg.create_mesh_XY(xRes=shape[0], yRes=shape[1])
    interpol_field = fg.interpolation_complex(field, mesh=mesh, fill_value=False)
    xMinMax = int(-shape[0] // 2 * xMinMax_frac[0]), int(shape[0] // 2 * xMinMax_frac[1])
    yMinMax = int(-shape[1] // 2 * yMinMax_frac[0]), int(shape[1] // 2 * yMinMax_frac[1])
    xyMesh_interpol = fg.create_mesh_XY(
        xRes=resolution[0], yRes=resolution[1],
        xMinMax=xMinMax, yMinMax=yMinMax)
    return interpol_field(*xyMesh_interpol), xyMesh_interpol


def one_plane_propagator(field, dz, stepsNumber_p, stepsNumber_m=None, n0=1, k0=1):
    """
    Double side propagation wrapper for fg.propagator_split_step_3D_linear
    :param field: 2D complex field
    :param dz: step along z
    :param stepsNumber_p: number of steps (forward, p - plus) [there is a chance it's confused with m direction)
    :param stepsNumber_m: number of steps (back, m - minus)
    :param n0: refractive index
    :param k0: wave number
    :return: 3D field
    """
    if stepsNumber_m is None:
        stepsNumber_m = stepsNumber_p
    fieldPropMinus = fg.propagator_split_step_3D_linear(field, dz=-dz, zSteps=stepsNumber_p, n0=n0, k0=k0)

    fieldPropPLus = fg.propagator_split_step_3D_linear(field, dz=dz, zSteps=stepsNumber_m, n0=n0, k0=k0)
    fieldPropTotal = np.concatenate((np.flip(fieldPropMinus, axis=2), fieldPropPLus[:, :, 1:]), axis=2)
    return fieldPropTotal


def main_field_processing(
        path,
        plotting=True,
        resolution_iterpol_center=(70, 70),
        xMinMax_frac_center=(1, 1),
        yMinMax_frac_center=(1, 1),
        resolution_interpol_working=(150, 150),
        xMinMax_frac_working=(1, 1),
        yMinMax_frac_working=(1, 1),
        resolution_crop=(120, 120),
        moments_init=None,
        moments_center=None,
):
    """
    This function:
     1) reading the field from matlab file
     2) converting it into numpy array
     3) normalizing
     4) finding the beam waste
     5) rescaling the field, using the interpolation, for faster next steps
     6) finding the beam center
     7) rescaling field to the scale we want for 3D calculations
     8) removing the tilt and shift

    Assumption
    ----------
    Beam waist finder only works with a uniform grid (dx = dy)

    :param path: file name
    :param plotting: if we want to see the plots and extra information
    :param resolution_iterpol_center: resolution for the beam center finder
    :param xMinMax_frac_center: rescale ration along X axis for the beam center
    :param yMinMax_frac_center: rescale ration along Y axis for the beam center
    :param resolution_interpol_working: resolution for the final field before the cropping
    :param xMinMax_frac_working: rescale ration along X axis for the beam center
    :param yMinMax_frac_working: rescale ration along X axis for the beam center
    :param resolution_crop: actual final resolution of the field
    :param moments_init: the moments for the LG spectrum
    :param moments_center: the moments for the beam center finder
    :return: 2D complex field
    """
    # beam width search work only with x_res==y_res

    if moments_init is None:
        moments_init = {'p': (0, 6), 'l': (-4, 4)}
    if moments_center is None:
        moments_center = {'p0': (0, 4), 'l0': (-4, 2)}

    # reading file
    field_init = read_field_2D_single(path)
    if plotting:
        plot_field(field_init)

    # normalization
    field_norm = normalization_field(field_init)
    if plotting:
        plot_field(field_norm)

    # creating mesh
    mesh_init = fg.create_mesh_XY(xRes=np.shape(field_norm)[0], yRes=np.shape(field_norm)[1])

    # finding beam waste
    width = find_beam_waist(field_norm, mesh=mesh_init)
    if plotting:
        print(f'Approximate beam waist: {width}')

    # rescaling field
    field_interpol, mesh_interpol = field_interpolation(
        field_norm, mesh=mesh_init, resolution=resolution_iterpol_center,
        xMinMax_frac=xMinMax_frac_center, yMinMax_frac=yMinMax_frac_center
    )
    if plotting:
        plot_field(field_interpol)

    # rescaling the beam width
    width = width / np.shape(field_norm)[0] * np.shape(field_interpol)[0]

    # plotting spec to select moments. .T because Danilo's code saving it like that
    if plotting:
        _ = cbs.LG_spectrum(field_interpol.T, **moments_init, mesh=mesh_interpol, plot=plotting, width=width, k0=1)

    # finding the beam center
    moments_init.update(moments_center)
    moments = moments_init
    x, y, eta, gamma = cbs.beamFullCenter(
        field_interpol, mesh_interpol,
        stepXY=(1, 1), stepEG=(3 / 180 * np.pi, 1 / 180 * np.pi),
        x=0, y=0, eta2=0., gamma=0.,
        **moments, threshold=1, width=width, k0=1, print_info=plotting
    )

    # rescaling field to the scale we want for 3D calculations
    field_interpol2, mesh_interpol2 = field_interpolation(
        field_norm, mesh=mesh_init, resolution=resolution_interpol_working,
        xMinMax_frac=xMinMax_frac_working, yMinMax_frac=yMinMax_frac_working
    )
    if plotting:
        plot_field(field_interpol2)

    # removing the tilt
    field_untilted = cbs.removeTilt(field_interpol2, mesh_interpol2, eta=-eta, gamma=gamma, k=1)
    if plotting:
        plot_field(field_untilted)

    # scaling the beam center
    shape = np.shape(field_untilted)
    x = int(x / np.shape(field_interpol)[0] * shape[0])
    y = int(y / np.shape(field_interpol)[1] * shape[1])

    # cropping the beam around the center
    field_cropped = field_untilted[
                    shape[0] // 2 - x - resolution_crop[0] // 2:shape[0] // 2 - x + resolution_crop[0] // 2,
                    shape[1] // 2 - y - resolution_crop[1] // 2:shape[1] // 2 - y + resolution_crop[1] // 2]
    if plotting:
        plot_field(field_cropped)

    # selecting the working field and field
    mesh = fg.create_mesh_XY(xRes=np.shape(field_cropped)[0], yRes=np.shape(field_cropped)[1])
    field = field_cropped
    print(f'field finished: {path[-20:]}')
    return field, mesh


if __name__ == '__main__':
    # test_hopf_turb_path = '3foil_turb_1.mat'
    # field2D, _ = main_field_processing(
    #     path=test_hopf_turb_path,
    #     plotting=True,
    #     resolution_iterpol_center=(60, 60),
    #     xMinMax_frac_center=(1, 1),
    #     yMinMax_frac_center=(1, 1),
    #     resolution_interpol_working=(150, 150),
    #     xMinMax_frac_working=(0.8, 0.8),
    #     yMinMax_frac_working=(0.8, 0.8),
    #     resolution_crop=(120, 120),
    #     moments_init={'p': (0, 6), 'l': (-4, 4)},
    #     moments_center={'p0': (0, 4), 'l0': (-4, 2)},
    # )
    # np.save('field2D_test.npy', field2D)
    field2D = np.load('field2D_test.npy')
    plot_field(field2D)

    # 2-directional propagation
    dz, steps_both = 5, 25
    field3D = one_plane_propagator(field2D, dz=dz, stepsNumber_p=steps_both, stepsNumber_m=None, n0=1, k0=1)
    # plot_field_3D_multi_planes(field3D, number=6, rows=3)
    plot_field(field3D[:, :, 0])
    plot_field(field3D[:, :, 25])
    plot_field(field3D[:, :, -1])

# x=0, y=0, eta=6.0*, gamma=2.0*, var=0.08804003185287904  70 70
