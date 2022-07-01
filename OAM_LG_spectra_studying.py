import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import numpy as np


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    xyMesh = fg.create_mesh_XY(xMinMax=3, yMinMax=3)
    beamInterpolated = interpolation_complex(beam, xArray, yArray)
    # beamInterpolated = beamInterpolated[0] + 1j * beamInterpolated[1]
    # field = beamInterpolated[0] + 1j * beamInterpolated[1]

    # pl.plot_2D(np.abs(beamInterpolated(*xyMesh)))
    # plt.show()
    rArray = np.linspace(0.01, 3, 25)
    aRArray = np.zeros(25, dtype=complex)
    for i in range(25):
        aRArray[i] = integral_number2_OAMcoefficients_FengLiPaper(beamInterpolated, rArray[i], 2)
    print(np.abs(aRArray))
