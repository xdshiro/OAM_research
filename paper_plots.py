import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
from my_functions.plotings import box_set_go
import numpy as np


def axis_set_go(fig, xyzMinMax=(-1, 1, -1, 1, -1, 1), mesh=None):
    fig.update_layout(font_size=24, font_family="Times New Roman", font_color='black',
                      legend_font_size=20,

                      scene=dict(
                          # annotations=[dict(x=[ 2], y=[-2], z=[-.5], ax=-2, ay=-2), dict(align='left'),
                          #              dict(align='left')],
                          xaxis_title=dict(text='', font=dict(size=30)),
                          yaxis_title=dict(text='', font=dict(size=30)),
                          zaxis_title=dict(text='', font=dict(size=30)),

                          aspectratio_x=2, aspectratio_y=2, aspectratio_z=2,

                          xaxis=dict(showticklabels=False, zeroline=False,
                                     showgrid=False, gridcolor='white',
                                     showbackground=False  # backgroundcolor='white',
                                     ),
                          yaxis=dict(showticklabels=False, zeroline=False,
                                     showgrid=False, gridcolor='white',
                                     showbackground=False
                                     ),
                          zaxis=dict(showticklabels=False, zeroline=False,
                                     showgrid=False, gridcolor='white',
                                     showbackground=False
                                     ), ),  # range=[-0.5, 0.5],
                      )
    xMin, xMax = xyzMinMax[0], xyzMinMax[1]
    yMin, yMax = xyzMinMax[2], xyzMinMax[3]
    zMin, zMax = xyzMinMax[4], xyzMinMax[5]
    lineX = np.array([[xMin, yMin, zMin], [xMax, yMin, zMin]])
    pl.plot_3D_dots_go(lineX, fig=fig, mode='lines', line={'width': 5, 'color': 'black'})
    labelX = np.array([[xMax, yMin, zMin]])
    pl.plot_3D_dots_go(labelX, fig=fig, mode='text', text=['x'])
    lineY = np.array([[xMin, yMin, zMin], [xMin, yMax, zMin]])
    pl.plot_3D_dots_go(lineY, fig=fig, mode='lines', line={'width': 5, 'color': 'black'})
    labelY = np.array([[xMin, yMax, zMin]])
    pl.plot_3D_dots_go(labelY, fig=fig, mode='text', text=['y'])
    lineZ = np.array([[xMin, yMin, zMin], [xMin, yMin, zMax]])
    pl.plot_3D_dots_go(lineZ, fig=fig, mode='lines', line={'width': 5, 'color': 'black'})
    labelZ = np.array([[xMin, yMin, zMax]])
    pl.plot_3D_dots_go(labelZ, fig=fig, mode='text', text=['z'])


if __name__ == '__main__':
    plot_10_examples = False
    if plot_10_examples:
        directoryName = f'C:\\Users\\Dima\\Box\\Knots Exp\\Experimental Data' \
                        f'\\dots\\trefoil\\Previous\\SR=0.95\\'
        fileName = '3foil_turb_35.npy'
        dotsDict = np.load(directoryName + fileName, allow_pickle=True).item()
        dots = sing.get_singularities(dotsDict)
        dots = np.array([dot for dot in dots if dot[2] < 50])
        newDots = sing.dots_filter(dots, checkValue=2, checkNumber=1)
        # filtering 2-4 dots clusters
        dots = sing.dots_filter(newDots, checkValue=4, checkNumber=3)
        print(dots)
        # pl.plot_scatter_3D(dots[:,0],dots[:,1],dots[:,2])
        fig = pl.plot_3D_dots_go(dots, marker={'size': 10, 'color': 'black', 'line': dict(width=185,
                                                                                              color='white')},
                                 )
        xMin, xMax = 1e10, 0
        yMin, yMax = 1e10, 0
        zMin, zMax = 1e10, 0
        for dot in dots:
            if dot[0] < xMin:
                xMin = dot[0]
            if dot[0] > xMax:
                xMax = dot[0]
            if dot[1] < yMin:
                yMin = dot[1]
            if dot[1] > yMax:
                yMax = dot[1]
            if dot[2] < zMin:
                zMin = dot[2]
            if dot[2] > zMax:
                zMax = dot[2]
        print(xMin, xMax, yMin, yMax, zMin, zMax)
        box_set_go(fig, xyzMinMax=(xMin, xMax, yMin, yMax, zMin, zMax), width=3)
        fig.show()

    scheme_plotting = False
    if scheme_plotting:
        reading_the_file = False
        if reading_the_file:
            fileName = f'C:\\Users\\Dima\\Box\\Knots Exp\\Experimental Data\\7-13-2022\\Field SR = 0.95\\3foil_turb_25.mat'
            field_experiment = fg.reading_file_mat(fileName=fileName, fieldToRead='U',
                                                   printV=False)
            from scipy.io import savemat

            mdic = {'field_exp': field_experiment}
            savemat('field_exp.mat', mdic)
            exit()

        # fieldAfterProp = fg.cut_filter(field_experiment, radiusPix=int(np.shape(field_experiment)[0] / 3.5), circle=True)
        # from scipy.io import savemat
        #
        # mdic = {'field_exp_cut': fieldAfterProp}
        # savemat('field_exp_cut.mat', mdic)
        # pl.plot_2D(np.abs(fieldAfterProp[:, :]))
        # pl.plot_2D(np.angle(fieldAfterProp[:, :]))

        # xMinMax = 2.0
        # yMinMax = 2.0
        # zMinMax = 0.5
        # zRes = 100
        # xRes = yRes = 200
        # xyzMesh = fg.create_mesh_XYZ(xMinMax, yMinMax, zMinMax, xRes, yRes, zRes, zMin=None)
        plot_field_3d = False
        if plot_field_3d:
            fileName = 'field_exp.mat'
            field_experiment = fg.reading_file_mat(fileName=fileName, fieldToRead='field_exp',
                                                   printV=False)
            fieldAfterProp = fg.one_plane_propagator(field_experiment, dz=10.5, stepsNumber=32, n0=1, k0=1)
            shape = np.shape(fieldAfterProp)
            xyCut = 20
            # pl.plot_2D(np.abs(fieldAfterProp[xyCut: shape[0] - xyCut,
            #                   xyCut: shape[1] - xyCut, 32]) ** 2)
            fieldAfterProp = fieldAfterProp[xyCut: shape[0] - xyCut, xyCut: shape[1] - xyCut, :]
            fieldAfterProp = fieldAfterProp / np.abs(fieldAfterProp).max()
            fig = pl.plot_3D_density(np.abs(fieldAfterProp), xMinMax=[-3, 3], yMinMax=[-3, 3], zMinMax=[-0.5, 0.5],
                                     opacityscale=[[0, 0], [0.2, 0.3], [0.3, 0.8], [1, 1]],
                                     # , [0.01, 0.3], [0.05, 0.6], [0.2, 0.6]
                                     surface_count=15, colorscale='jet', show=False)
            box_set_go(fig, xyzMinMax=(-3, 3, -3, 3, -0.5, 0.5), width=2)
            # fig.update_layout(font_size=18, font_family="Times New Roman", font_color='black',
            #                   legend_font_size=20,
            #                   scene=dict(
            #                       xaxis_title=dict(text='x', font=dict(size=30)),
            #                       yaxis_title=dict(text='y', font=dict(size=30)),
            #                       zaxis_title=dict(text='z', font=dict(size=30)),
            #
            #                       aspectratio_x=2, aspectratio_y=2, aspectratio_z=2,
            #                       xaxis=dict(range=[-2, 2], showticklabels=False, nticks=7),
            #                       # # range=[-xMinMax, xMinMax]
            #                       yaxis=dict(range=[-2, 2], showticklabels=False, nticks=7),  # range=[-2, 2]
            #                       zaxis=dict(nticks=7, showticklabels=False), ),  # range=[-0.5, 0.5],
            #                   )

            fig.show()
            exit()
        plot_one_plane_dots = False
        if plot_one_plane_dots:
            fileName = 'field_exp.mat'
            field_experiment = fg.reading_file_mat(fileName=fileName, fieldToRead='field_exp',
                                                   printV=False)
            fieldAfterProp = fg.one_plane_propagator(field_experiment, dz=7.5, stepsNumber=32, n0=1, k0=1)
            shape = np.shape(fieldAfterProp)
            fieldAfterProp = fg.cut_filter(fieldAfterProp, radiusPix=int(np.shape(fieldAfterProp)[0] // 3.5),
                                           circle=True)
            xyCut = 55
            # pl.plot_2D(np.abs(fieldAfterProp[xyCut: shape[0] - xyCut,
            #                   xyCut: shape[1] - xyCut, 32]) ** 2)
            fieldAfterProp = fieldAfterProp[xyCut + 10: shape[0] - xyCut, xyCut: shape[1] - xyCut, :]
            # pl.plot_2D(np.abs(fieldAfterProp[:, :, 32]))
            # exit()
            fieldAfterProp = fieldAfterProp[20:85, 20:90, :]
            dots = sing.get_singularities(np.angle(fieldAfterProp), bigSingularity=False, axesAll=False)
            fig = pl.plot_3D_dots_go(dots, marker={'size': 10, 'color': 'black', 'line': dict(width=3,
                                                                                                  color='grey')}, show=False)
            box_set_go(fig, xyzMinMax=(0, 65, 0, 70, 0, shape[2] - 4), width=2)
            # axis_set_go(fig, xyzMinMax=(-1, 1, -1, 1, -1, 1))
            fig.show()
            exit()
            # fig = pl.plot_3D_density(np.abs(fieldAfterProp), xMinMax=[-3, 3], yMinMax=[-3, 3], zMinMax=[-0.5, 0.5],
            #                          opacityscale=[[0, 0], [0.2, 0.3], [0.3, 0.8], [1, 1]],
            #                          # , [0.01, 0.3], [0.05, 0.6], [0.2, 0.6]
            #                          surface_count=15, colorscale='jet', show=False)
            # fig.show()

        plot_t_plane_dots = True
        if plot_t_plane_dots:
            fileName = 'field_exp.mat'
            field_experiment = fg.reading_file_mat(fileName=fileName, fieldToRead='field_exp',
                                                   printV=False)
            fieldAfterProp = fg.one_plane_propagator(field_experiment, dz=7.5, stepsNumber=32, n0=1, k0=1)
            shape = np.shape(fieldAfterProp)
            fieldAfterPropIn = fg.cut_filter(fieldAfterProp, radiusPix=int(np.shape(fieldAfterProp)[0] // 3.5),
                                           circle=True)
            xyCut = 55
            # pl.plot_2D(np.abs(fieldAfterProp[xyCut: shape[0] - xyCut,
            #                   xyCut: shape[1] - xyCut, 32]) ** 2)
            width = 3
            fieldAfterProp = fieldAfterPropIn[xyCut + 10: shape[0] - xyCut, xyCut: shape[1] - xyCut, :]
            # fieldAfterProp = fieldAfterProp[20:85, 20:90, :]
            # pl.plot_2D(np.abs(fieldAfterProp[:, :, 32]))
            dots = sing.get_singularities(np.angle(fieldAfterProp), bigSingularity=False, axesAll=True)
            fig = pl.plot_3D_dots_go(dots, marker={'size': 10, 'color': 'skyblue', 'line': dict(width=width,
                                                                                             color='black')})
            newDots = sing.dots_filter(dots, checkValue=2, checkNumber=1)
            # filtering 2-4 dots clusters
            newDots2 = sing.dots_filter(newDots, checkValue=4, checkNumber=3)
            pl.plot_3D_dots_go(newDots2, marker={'size': 10, 'color': 'black', 'line': dict(width=width,
                                                                                            color='grey')}, fig=fig)
            box_set_go(fig, xyzMinMax=(0, 105, 0, 105, 0, shape[2] - 4), width=2)
            fig.show()
            fieldAfterProp = fieldAfterProp[20:85, 20:90, :]
            dots = sing.get_singularities(np.angle(fieldAfterProp), bigSingularity=False, axesAll=True)
            newDots = sing.dots_filter(dots, checkValue=2, checkNumber=1)
            # filtering 2-4 dots clusters
            newDots2 = sing.dots_filter(newDots, checkValue=4, checkNumber=3)
            fig = pl.plot_3D_dots_go(newDots2, marker={'size': 10, 'color': 'skyblue', 'line': dict(width=width,
                                                                                                 color='black')})
            newDots3 = sing.dots_dens_reduction(newDots2, checkValue=2.5, checkNumber=3)
            pl.plot_3D_dots_go(newDots3, marker={'size': 10, 'color': 'black', 'line': dict(width=width,
                                                                                            color='grey')}, fig=fig)
            box_set_go(fig, xyzMinMax=(0, 65, 0, 70, 0, shape[2] - 4), width=2)
            fig.show()
            fig = pl.plot_3D_dots_go(newDots3, marker={'size': 10, 'color': 'grey', 'line': dict(width=width,
                                                                                                 color='black')})
            fig = pl.plot_3D_dots_go(newDots3, marker={'size': 10, 'color': 'black', 'line': dict(width=width,
                                                                                                  color='grey')},
                                     fig=fig)
            box_set_go(fig, xyzMinMax=(0, 65, 0, 70, 0, shape[2] - 4), width=2)
            fig.show()
            exit()
            # fig = pl.plot_3D_density(np.abs(fieldAfterProp), xMinMax=[-3, 3], yMinMax=[-3, 3], zMinMax=[-0.5, 0.5],
            #                          opacityscale=[[0, 0], [0.2, 0.3], [0.3, 0.8], [1, 1]],
            #                          # , [0.01, 0.3], [0.05, 0.6], [0.2, 0.6]
            #                          surface_count=15, colorscale='jet', show=False)
            # fig.show()


        plot_pyknotid = False
        if plot_pyknotid:
            fileName = 'field_exp.mat'
            field_experiment = fg.reading_file_mat(fileName=fileName, fieldToRead='field_exp',
                                                   printV=False)
            fieldAfterProp = fg.one_plane_propagator(field_experiment, dz=7.5, stepsNumber=32, n0=1, k0=1)
            shape = np.shape(fieldAfterProp)
            fieldAfterProp = fg.cut_filter(fieldAfterProp, radiusPix=int(np.shape(fieldAfterProp)[0] // 3.5),
                                           circle=True)
            xyCut = 55
            # pl.plot_2D(np.abs(fieldAfterProp[xyCut: shape[0] - xyCut,
            #                   xyCut: shape[1] - xyCut, 32]) ** 2)
            width = 3
            fieldAfterProp = fieldAfterProp[xyCut + 10: shape[0] - xyCut, xyCut: shape[1] - xyCut, :]
            # pl.plot_2D(np.abs(fieldAfterProp[:, :, 32]))
            dots = sing.get_singularities(np.angle(fieldAfterProp), bigSingularity=False, axesAll=True)
            dotsKnot = sing.knot_sequence_from_dots(dots, checkValue1=2, checkNumber1=1,
                                               checkValue2=4, checkNumber2=3,
                                               checkValue3=2.8, checkNumber3=3)
            # trefoil = sp.Knot(dotsKnot, add_closure=False)
            trefoil = sing.knot_build_pyknotid(dotsKnot)
            trefoil.plot(tube_radius=2.)
            from vispy import app

            app.run()
            trefoil.interpolate(250, quiet=True, per=False)
            trefoilSmooth = trefoil.points[2: -1]
            import pyknotid.spacecurves as sp
            trefoil = sp.Knot(trefoilSmooth, add_closure=False)
            trefoil.plot(tube_radius=2.)
            app.run()
        # pl.plot_2D(np.abs(field_experiment))
        # pl.plot_2D(np.angle(field_experiment))

        exit()
        # dotsExp = np.load('C:\\Users\\Dima\\Box\\Knots Exp\\Experimental Data\\dots\\trefoil\\'
        #                    'Field 3foil noturb\\3foil_noturb_1.npy',  # 25
        #                    allow_pickle=True).item()
        dotsExp = np.load('C:\\Users\\Dima\\Box\\Knots Exp\\Experimental Data\\dots\\trefoil\\Previous\\'
                          'Field SR = 0.95\\3foil_turb_25.npy',  # 25
                          allow_pickle=True).item()
        dots = sing.get_singularities(dotsExp)
        # dotsKnot = sing.knot_sequence_from_dots(dots, checkValue1=2, checkNumber1=1,
        #                                                                                 checkValue2=4, checkNumber2=3,
        #                                                                                 checkValue3=3, checkNumber3=3)
        # pl.plot_scatter_3D(dotsKnot[:, 0], dotsKnot[:, 1], dotsKnot[:, 2], size=100)
        # dotsKnot = sing.knot_sequence_from_dots(dots, checkValue1=2, checkNumber1=1,
        #                                    checkValue2=4, checkNumber2=3,
        #                                    checkValue3=3, checkNumber3=3)
        newDots = sing.dots_filter(dots, checkValue=2, checkNumber=1)
        # filtering 2-4 dots clusters
        newDots2 = sing.dots_filter(newDots, checkValue=4, checkNumber=3)
        newDots3 = sing.dots_dens_reduction(newDots2, checkValue=4, checkNumber=3)
        dotsKnot = newDots3
        # dotsKnot = sing.knot_sequence_from_dots(dots, checkValue1=2, checkNumber1=1,
        #                                         checkValue2=4, checkNumber2=3,
        #                                         checkValue3=0, checkNumber3=3)
        fig = pl.plot_3D_dots_go(dotsKnot, marker={'size': 10, 'color': 'black'})

        exit()
        sing.plot_knot_pyknotid(dotsKnot, interpolation=300, tube_radius=2.0, per=False, add_closure=False,
                                tube_points=14, fov=0, flip=(False, False, True))

    crossections_on_knot = False
    if crossections_on_knot:
        xMinMax = 2.1
        yMinMax = 2.1
        zMinMax = 0.5
        zRes = 12
        xRes = yRes = 24
        xyzMesh = fg.create_mesh_XYZ(xMinMax, yMinMax, zMinMax, xRes, yRes, zRes, zMin=None)
        modesTrefoil = [(0, 0), (0, 1), (0, 2), (0, 3), (3, 0)]
        coeffTrefoil = [1.29, -3.95, 7.49, -3.28, -3.98]
        beam = bp.LG_combination(*xyzMesh, coefficients=coeffTrefoil, modes=modesTrefoil)
        dots = sing.get_singularities(np.angle(beam), bigSingularity=False, axesAll=False)
        dots = fg.dots3D_rescale(dots, xyzMesh)
        beamToPlot = np.abs(beam) / np.abs(beam).max()
        width = 3
        fig = pl.plot_3D_dots_go(dots, marker={'size': 10, 'color': 'black', 'line': dict(width=width,
                                                                                             color='grey')})
        xyz = fg.arrays_from_mesh(xyzMesh)
        z1, z2 = zRes / 1.17, zRes / 2
        zPlane1 = xyz[2][int(z1)] * np.ones((xRes, yRes))
        zPlane2 = xyz[2][int(z2)] * np.ones((xRes, yRes))
        extraDotsRed = np.array([list(dot) for dot in dots if (dot[2] == zPlane1[0, 0] or dot[2] == zPlane2[0, 0])])
        pl.plot_plane_go(zPlane1, xyzMesh, fig=fig)
        pl.plot_plane_go(zPlane2, xyzMesh, fig=fig)
        pl.plot_3D_dots_go(extraDotsRed, marker={'size': 16, 'color': 'red', 'line': dict(width=width,
                                                                                             color='black')}, fig=fig)
        # pl.plot_3D_density(beamToPlot, mesh=xyzMesh, fig=fig,
        #                          colorscale='jet',
        #                          opacityscale=[[0, 0], [0.2, 0], [0.3, 0.5], [1, 1]], opacity=0.4, surface_count=20,
        #                          show=False)
        mycolorscale = [[0, '#aa9ce2'],
                        [1, '#aa9ce2']]
        print(extraDotsRed)
        # xMult = 1.0
        # yMult = 1.0
        # zMult = 1.0
        # lines = np.array([[-xMinMax  ,-yMinMax,-zMinMax], [xMinMax * xMult,-yMinMax,-zMinMax]])
        # pl.plot_3D_dots_go(lines, fig=fig, mode='lines', line={'width': 5, 'color': 'black'})
        # pl.plot_3D_dots_go(lines, fig=fig, mode='text', text=['', 'x'])
        # lines = np.array([[-xMinMax, -yMinMax, -zMinMax], [-xMinMax, yMinMax * yMult, -zMinMax]])
        # pl.plot_3D_dots_go(lines, fig=fig, mode='lines', line={'width': 5, 'color': 'black'})
        # pl.plot_3D_dots_go(lines, fig=fig, mode='text', text=['', 'y'])
        # lines = np.array([[-xMinMax, -yMinMax, -zMinMax], [-xMinMax, -yMinMax, zMinMax * zMult]])
        # pl.plot_3D_dots_go(lines, fig=fig, mode='lines', line={'width': 5, 'color': 'black'}, text=['a', 'b'])
        # lines = np.array([[-xMinMax * 0.9, -yMinMax, zMinMax * 0.9]])
        # pl.plot_3D_dots_go(lines, fig=fig, mode='text', text=['z'])

        # fig.update_layout(
        #     scene=dict(
        #         aspectratio_x=2, aspectratio_y=2, aspectratio_z=2,
        #         xaxis=dict(range=[-xMinMax, xMinMax], nticks=7),
        #         yaxis=dict(range=[-yMinMax, yMinMax]),
        #         zaxis=dict(range=[-zMinMax, zMinMax]), ),
        #     # width=700,
        #     # margin=dict(r=0, l=10, b=70, t=10)
        # )
        # ann = [dict(x=x, y=y, z=z, text='F') for x, y, z in zip(anx, any, anz)]
        box_set_go(fig, xyzMinMax=(-xMinMax, xMinMax, -yMinMax, yMinMax, -zMinMax, zMinMax), width=2)
        # fig.update_layout(font_size=24, font_family="Times New Roman", font_color='black',
        #                   legend_font_size=20,
        #
        #                   scene=dict(
        #                       # annotations=[dict(x=[ 2], y=[-2], z=[-.5], ax=-2, ay=-2), dict(align='left'),
        #                       #              dict(align='left')],
        #                       xaxis_title=dict(text='x', font=dict(size=30)),
        #                       yaxis_title=dict(text='y', font=dict(size=30)),
        #                       zaxis_title=dict(text='z', font=dict(size=30)),
        #
        #                       aspectratio_x=2, aspectratio_y=2, aspectratio_z=2,
        #
        #                       xaxis=dict(range=[-xMinMax, xMinMax], nticks=5,
        #                                  showticklabels=False, gridcolor='white',
        #                                  showgrid=True, zeroline=False, showbackground=True,  # backgroundcolor='white',
        #                                  ),
        #                       yaxis=dict(range=[-yMinMax, yMinMax], nticks=5
        #                                  , showticklabels=False, gridcolor='white',
        #                                  showgrid=True, zeroline=False, showbackground=True),  # range=[-2, 2]
        #                       zaxis=dict(range=[-zMinMax, zMinMax], nticks=5
        #                                  , showticklabels=False, gridcolor='white',
        #                                  showgrid=True, zeroline=False, showbackground=True), ),  # range=[-0.5, 0.5],
        #                   )
        fig.show()
        exit()

    plot_trefoil_w = True
    if plot_trefoil_w:
        xMinMax = 3.8
        yMinMax = 3.8
        zMinMax = 2
        zRes = 3
        xRes = yRes = 512
        xyzMesh = fg.create_mesh_XYZ(xMinMax, yMinMax, zMinMax, xRes, yRes, zRes, zMin=None) # , xMin= -3
        trefoil = bp.trefoil(*xyzMesh, w=1.5, width=1.0, k0=1.0)
        # trefoil = bp.hopf(*xyzMesh, w=1.376, width=1.0, k0=1.0)
        from scipy.io import savemat
        #
        mdic = {"field": trefoil}
        savemat('trefoil15.mat', mdic)
        pl.plot_2D(np.angle(trefoil[:, :, zRes//2]), show=True)
        pl.plot_2D(np.abs(trefoil[:, :, zRes//2]), show=True)
        exit()
        dots = sing.get_singularities(np.angle(trefoil), axesAll=True)
        # pl.plot_scatter_3D(dots[:, 0], dots[:, 1], dots[:, 2])
        fig = pl.plot_3D_dots_go(dots, show=False, marker={'size': 10, 'color': 'black', 'line': dict(width=3,
                                                                                             color='grey')})
        box_set_go(fig=fig, xyzMinMax=(0, xRes, 0, yRes, 0, zRes))
        fig.show()