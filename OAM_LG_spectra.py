import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np
import matplotlib.pyplot as plt

# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    def beam_comination(x, y):
        return bp.LG_combination(x, y, coefficients=[1, 1, 1], modes=[(1, 0), (2, 0), (3, 0)])


    def beam_LG(x, y):
        return bp.LG_simple(x, y, l=2, width=1.3)


    file = f'./DATA/OAM_LG_spectra_studying/v1_field_SR_6.000000e-01_num_2.mat'
    beamTurb = fg.reading_file_mat(file, fieldToRead='Efield')[:, :, 0]

    xyMesh = fg.create_mesh_XY(xMax=10, yMax=10, xRes=np.shape(beamTurb)[0], yRes=np.shape(beamTurb)[1])
    xyArrays = fg.arrays_from_mesh(xyMesh)
    # beam = bp.LG_combination(*xyMesh, [1, 1, 1], [(1, 0), (2, 0), (3, 0)])
    beamTurbInter = fg.interpolation_complex(beamTurb, *xyArrays)
    # pl.plot_2D(np.abs(beamTurbInter(*xyMesh)))
    # pl.plot_2D(np.angle(beamTurbInter(*xyMesh)))
    # plt.show()
    # pl.plot_2D(np.abs(beam_LG(*xyMesh)))
    # pl.plot_2D(np.angle(beam_LG(*xyMesh)))
    # plt.show()
    # beam = beam_comination(*xyMesh)
    # pl.plot_2D(np.abs(beam))
    # plt.show()
    # exit()
    # C:\WORK\CODES\OAM_research\DATA\OAL_LG_spectra_studying
    # v1_SR_7.000000e-01_num_1
    rArray = np.linspace(0.01, 10, 25)
    aRArray = np.zeros(25, dtype=complex)
    from scipy import integrate

    for l in [-1, -2, -3]:
        for i in range(25):
            aRArray[i] = sing.integral_number2_OAMcoefficients_FengLiPaper(beamTurbInter, rArray[i], l)
        print(integrate.simps(np.abs(aRArray) ** 2 * rArray, rArray))
        plt.plot(rArray, np.imag(aRArray), label=f'l={abs(l)}')
    plt.legend()
    plt.title(f'Im(a_l)')
    plt.show()
