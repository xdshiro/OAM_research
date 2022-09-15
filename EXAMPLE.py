import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np

if __name__ == '__main__':
    plot_trefoil = True
    if plot_trefoil:
        xyzMesh = fg.create_mesh_XYZ(2.1, 2.1, 0.7, 50, 50, 50, zMin=None)
        beam = bp.LG_combination(*xyzMesh,
                                 coefficients=[1.71, -5.66, 6.38, -2.30, -4.36],
                                 modes=[(0, 0), (0, 1), (0, 2), (0, 3), (3, 0)],
                                 width=[1, 1, 1, 1, 1])
        dots = sing.get_singularities(np.angle(beam))
        # pl.plot_2D(np.abs(beam[:, :, zRes//2]))
        # pl.plot_scatter_3D(dots[:, 0], dots[:, 1], dots[:, 2])
        fig = pl.plot_3D_dots_go(dots)
        pl.box_set_go(fig, mesh=None, autoDots=dots, perBox=0.05)
        fig.show()
