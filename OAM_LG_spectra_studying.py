import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np
import matplotlib.pyplot as plt



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    xyMesh = fg.create_mesh_XY(xMax=4, yMax=4)
    xyArrays = fg.arrays_from_mesh(xyMesh)
    def beam_comination(x, y):
        return bp.LG_combination(x, y, coefficients=[1, 1, 1], modes=[(1, 0), (2, 0), (3, 0)])
    # beam = bp.LG_combination(*xyMesh, [1, 1, 1], [(1, 0), (2, 0), (3, 0)])
    # beamInterpolated = fg.interpolation_complex(beam, *xyArrays)
    # beam = beam_comination(*xyMesh)
    # pl.plot_2D(np.abs(beam))
    # plt.show()
    # exit()
    # C:\WORK\CODES\OAM_research\DATA\OAL_LG_spectra_studying
    # v1_SR_7.000000e-01_num_1
    rArray = np.linspace(0.01, 4, 25)
    aRArray = np.zeros(25, dtype=complex)
    for l in [1, 2, 3]:
        for i in range(25):
            aRArray[i] = fg.integral_number2_OAMcoefficients_FengLiPaper(beam_comination, rArray[i], l)
        plt.plot(rArray, np.abs(aRArray))
    plt.show()