import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np

if __name__ == '__main__':
    xMinMax = 4
    yMinMax = 4
    zMinMax = 1
    zRes = 250
    xRes = yRes = 250
    xyzMesh = fg.create_mesh_XYZ(xMinMax, yMinMax, zMinMax, xRes, yRes, zRes, zMin=None)
    # modes = [(0, 0), (0, 1), (0, 2), (1, 0), (2, 0), (2, 1), (2, 2)]
    # modes = [(0, 0), (0, 1), (1, 0), (1, 1)]
    # lpArray = [0, 0, 0, 0]
    # diapason = [1, 1, 1, 1]
    modes = [(0, 0), (0, 1), (1, 0), (1, 1), (-1, 0), (-1, 1)]
    lpArray = [0, 0, 0, 0, 0, 0]
    diapason = [1, 1, 1, 1, 1, 1]
    diapasonComplex = lpArray
    while True:
        coeff = fg.random_list(lpArray, diapason, diapason_complex=diapasonComplex)
        print(coeff)
        beam = bp.LG_combination(*xyzMesh, coefficients=coeff, modes=modes)
        pl.plot_2D(np.abs(beam[:, :, zRes // 2]))
        sing.plot_knot_dots(beam, show=True, axesAll=False)

    exit()
    coeff = [1.715, -5.662, 6.381 / 16 * 17, -2.305 * 3, -4.356]
    phase = [0, 0, np.pi / 16 * 0, 0, 0]
    coeff = [a * np.exp(1j * p) for a, p in zip(coeff, phase)]
    beam = bp.LG_combination(*xyzMesh, coefficients=coeff, modes=((0, 0), (0, 1), (0, 2), (0, 3), (3, 0)))
    sing.plot_knot_dots(beam, show=True)
