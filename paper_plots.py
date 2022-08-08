import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np

if __name__ == '__main__':
    crossections_on_knot = True
    if crossections_on_knot:
        xMinMax = 2.1
        yMinMax = 2.1
        zMinMax = 0.5
        zRes = 100
        xRes = yRes = 200
        xyzMesh = fg.create_mesh_XYZ(xMinMax, yMinMax, zMinMax, xRes, yRes, zRes, zMin=None)
        modesTrefoil = [(0, 0), (0, 1), (0, 2), (0, 3), (3, 0)]
        coeffTrefoil = [1.29, -3.95, 7.49, -3.28, -3.98]
        beam = bp.LG_combination(*xyzMesh, coefficients=coeffTrefoil, modes=modesTrefoil)
        dots = sing.get_singularities(np.angle(beam), bigSingularity=False, axesAll=False)
        dots = fg.dots3D_rescale(dots, xyzMesh)
        beamToPlot = np.abs(beam) / np.abs(beam).max()
        fig = pl.plot_3D_dots_go(dots, marker={'size': 8, 'color': 'black'})
        xyz = fg.arrays_from_mesh(xyzMesh)
        z1, z2 = zRes / 1.17, zRes / 2
        zPlane1 = xyz[2][int(z1)] * np.ones((xRes, yRes))
        zPlane2 = xyz[2][int(z2)] * np.ones((xRes, yRes))
        extraDotsRed = np.array([list(dot) for dot in dots if (dot[2] == zPlane1[0, 0] or dot[2] == zPlane2[0, 0])])
        pl.plot_plane_go(zPlane1, xyzMesh, fig=fig)
        pl.plot_plane_go(zPlane2, xyzMesh, fig=fig)
        fig = pl.plot_3D_dots_go(extraDotsRed, marker={'size': 14, 'color': 'red'}, fig=fig)
        # pl.plot_3D_density(beamToPlot, mesh=xyzMesh, fig=fig,
        #                          colorscale='jet',
        #                          opacityscale=[[0, 0], [0.2, 0], [0.3, 0.5], [1, 1]], opacity=0.4, surface_count=20,
        #                          show=False)
        mycolorscale = [[0, '#aa9ce2'],
                        [1, '#aa9ce2']]

        fig.update_layout(
            scene=dict(
                aspectratio_x=2, aspectratio_y=2, aspectratio_z=2,
                xaxis=dict(range=[-xMinMax, xMinMax]),
                yaxis=dict(range=[-yMinMax, yMinMax]),
                zaxis=dict(range=[-zMinMax, zMinMax]), ),
            # width=700,
            # margin=dict(r=0, l=10, b=70, t=10)
        )
        fig.show()
        exit()
