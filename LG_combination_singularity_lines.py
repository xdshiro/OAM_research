import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np

if __name__ == '__main__':
    xMinMax = 3
    yMinMax = 3
    zMinMax = 1
    zRes = 80
    xRes = yRes = 80
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
        sing.plot_knot_dots(beam, show=True, axesAll=True)

    coeff = [(-0.7465130995545659+0j), (0.9634448817109214+0j), (-0.4191932274252561+0j), (0.12832580982839925+0j), (-0.5350233092389178+0j), (-0.7291809617196814+0j)]
    beam = bp.LG_combination(*xyzMesh, coefficients=coeff, modes=modes)
    sing.plot_knot_dots(beam, show=True, axesAll=True)
