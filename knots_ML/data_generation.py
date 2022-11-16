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


def find_beam_waist(field):
    """
    wrapper for the beam waste finder. More details in knots_ML.center_beam_search
    """
    shape = np.shape(field)
    xyMesh_init = fg.create_mesh_XY(xRes=shape[0], yRes=shape[1])
    width = cbs.find_width(field, mesh=xyMesh_init, width=shape[1] // 8, widthStep=1, print_steps=False)
    return width



if __name__ == '__main__':
    test_hopf_turb_path = ('C:\\Users\\Cmex-\\Box\\Knots Exp\\'
                           'Experimental Data\\7-13-2022\\Field SR = 0.85\\3foil_turb_1.mat')
    field_init = read_field_2D_single(test_hopf_turb_path)
    plot_field(field_init)
    field_norm = normalization_field(field_init)
    plot_field(field_norm)
    width = find_beam_waist(field_norm)
    print(f'Approximate beam waist: {width}')
