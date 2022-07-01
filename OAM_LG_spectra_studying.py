import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np
import matplotlib.pyplot as plt



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # xyMesh = fg.create_mesh_XY(xMax=4, yMax=3, xRes=4, yRes=3)
    # print(xyMesh)
    # exit()
    xyMesh = fg.create_mesh_XY(xMax=4, yMax=3, xRes=17, yRes=19)
    # print(len(np.shape(xyMesh[1])))
    # print(xyMesh[1][0, :, 1])
    beam = bp.trefoil(*xyMesh, w=1.4)
    beamInterpolated = fg.interpolation_complex(beam, *fg.arrays_from_mesh(xyMesh))
    pl.plot_2D(np.abs(beam), axis_equal=True)
    pl.plot_2D(np.abs(beamInterpolated(*xyMesh)), axis_equal=True)
    plt.show()
    exit()
    plt.show()
    # xyMesh
    exit()
    # print(np.shape(xyMesh[0][:, 0]), xyMesh[0][1], xyMesh[1][1])
    # exit()
    xyzMesh = fg.create_mesh_XYZ(xMax=4, yMax=3, zMax=1)
    print(xyzMesh[0][1][2])#, xyzMesh[1][1], xyzMesh[2][1])
    exit()

    # beamInterpolated = beamInterpolated[0] + 1j * beamInterpolated[1]
    # field = beamInterpolated[0] + 1j * beamInterpolated[1]

    # pl.plot_2D(np.abs(beamInterpolated(*xyMesh)))
    # plt.show()
    rArray = np.linspace(0.01, 3, 25)
    aRArray = np.zeros(25, dtype=complex)
    for i in range(25):
        aRArray[i] = integral_number2_OAMcoefficients_FengLiPaper(beamInterpolated, rArray[i], 2)
    print(np.abs(aRArray))
