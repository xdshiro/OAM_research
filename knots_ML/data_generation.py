"""

"""
import my_functions.functions_general as fg
# import importlib
# importlib.reload(fg)
import my_functions.singularities as sing
import my_functions.beams_and_pulses as bp
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


def one_plane_propagator(fieldPlane, dz, stepsNumber_p, stepsNumber_m=None, n0=1, k0=1):  # , shapeWrong=False
    if stepsNumber_m is None:
        stepsNumber_m = stepsNumber_p
    fieldPropMinus = fg.propagator_split_step_3D_linear(fieldPlane, dz=-dz, zSteps=stepsNumber_p, n0=n0, k0=k0)

    fieldPropPLus = fg.propagator_split_step_3D_linear(fieldPlane, dz=dz, zSteps=stepsNumber_m, n0=n0, k0=k0)
    fieldPropTotal = np.concatenate((np.flip(fieldPropMinus, axis=2), fieldPropPLus[:, :, 1:-1]), axis=2)
    return fieldPropTotal


def main_field_processing(
        path,
        plotting=True,
        resolution_iterpol=(70, 70),
        resolution_crop=(50, 50),
        moments_init=None,
        moments_center=None,
        ):
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
        field_norm, mesh=mesh_init, resolution=resolution_iterpol,
        xMinMax_frac=(1., 1.), yMinMax_frac=(1., 1.)
    )
    if plotting:
        plot_field(field_interpol)

    # rescaling the beam width
    width = width / np.shape(field_norm)[0] * np.shape(field_interpol)[0]

    # plotting spec to select moments. .T because Danilo's code saving it like that
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

    # removing the tilt
    field_untilted = cbs.removeTilt(field_interpol, mesh_interpol, eta=-eta, gamma=gamma, k=1)
    if plotting:
        plot_field(field_untilted)

    # cropping the beam around the center
    shape = np.shape(field_untilted)
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
    # test_hopf_turb_path = ('C:\\Users\\Cmex-\\Box\\Knots Exp\\'
    #                        'Experimental Data\\7-13-2022\\Field SR = 0.85\\3foil_turb_1.mat')
    # field2D, _ = main_field_processing(
    #     path=test_hopf_turb_path,
    #     plotting=True,
    #     resolution_iterpol=(70, 70),
    #     resolution_crop=(50, 50),
    #     moments_init={'p': (0, 6), 'l': (-4, 4)},
    #     moments_center={'p0': (0, 4), 'l0': (-4, 2)},
    # )
    # np.save('field2D_test.npy', field2D)
    field2D = np.load('field2D_test.npy')
    plot_field(field2D)
