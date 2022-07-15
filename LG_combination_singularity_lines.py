import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np

if __name__ == '__main__':
    xMinMax = 3
    yMinMax = 3
    zMinMax = 1
    zRes = 50
    xRes = yRes = 50
    xyzMesh = fg.create_mesh_XYZ(xMinMax, yMinMax, zMinMax, xRes, yRes, zRes, zMin=None)
    # modes = [(0, 0), (0, 1), (0, 2), (1, 0), (2, 0), (2, 1), (2, 2)]
    # modes = [(0, 0), (0, 1), (1, 0), (1, 1)]
    # lpArray = [0, 0, 0, 0]
    # diapason = [1, 1, 1, 1]
    modes = [(0, 0), (0, 1), (1, 0), (1, 1), (-1, 0), (-1, 1)]
    lpArray = [0, 0, 0, 0, 0, 0]
    diapason = [1, 1, 1, 1, 1, 1]
    diapasonComplex = lpArray
    while False:
        coeff = fg.random_list(lpArray, diapason, diapason_complex=diapasonComplex)
        print(coeff)
        beam = bp.LG_combination(*xyzMesh, coefficients=coeff, modes=modes)
        pl.plot_2D(np.abs(beam[:, :, zRes // 2]))
        sing.plot_knot_dots(beam, show=True, axesAll=True)

    coeff = [(-0.7465130995545659 + 0j), (0.9634448817109214 + 0j), (-0.4191932274252561 + 0j),
             (0.12832580982839925 + 0j), (-0.5350233092389178 + 0j), (-0.7291809617196814 + 0j)]
    beam = bp.LG_combination(*xyzMesh, coefficients=coeff, modes=modes)
    # sing.plot_knot_dots(beam, show=True, axesAll=True)
    beamToPlot = np.abs(beam) / np.abs(beam).max()
    import plotly.graph_objects as go

    # function to create fast dotsOnly, we also need scaling mesh
    dotsFull, dotsOnly = sing.cut_non_oam(np.angle(beam),
                                          bigSingularity=False, axesAll=False)
    dots = np.array([list(dots) for (dots, OAM) in dotsOnly.items()])

    # dotsMinus = np.array([list(dots) for (dots, OAM) in dotsOnly.items() if OAM == -1])
    fig = pl.plot_3D_density(np.abs(beam) / np.abs(beam).max(), colorscale='jet',
                             xMinMax=[0, xRes], yMinMax=[0, yRes], zMinMax=[0, zRes],
                             opacityscale=[[0, 0], [0.1, 0], [1, 1]], opacity=0.6, surface_count=20, show=False)
    fig.add_trace(go.Scatter3d(x=dots[:, 0], y=dots[:, 1], z=dots[:, 2],
                               mode='markers', marker = dict( color='rgb(0,0,0)')))
    # fig.add_trace(go.Scatter3d(x=[1], y=[0], z=[0],
    #                                    mode='markers'))

    # plt.show()
    fig.show()
    exit()
    # pl.plot_3D_density(beamToPlot, colorscale='jet',
    #                    opacity=0.0, surface_count=20, show=True)
