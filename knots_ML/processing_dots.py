import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np

def find_the_closest_dot(array):
    pass
if __name__ == '__main__':
    plot_trefoil = True
    if plot_trefoil:
        xyzMesh = fg.create_mesh_XYZ(2.1, 2.1, 0.7, 40, 40, 40, zMin=None)
        beam = bp.LG_combination(*xyzMesh,
                                 coefficients=[1.71, -5.66, 6.38, -2.30, -4.36],
                                 modes=[(0, 0), (0, 1), (0, 2), (0, 3), (3, 0)],
                                 width=[1, 1, 1, 1, 1])
        pl.plot_2D(np.abs(beam[:,:,40//2]))
        sing.Jz_calc_no_conj(beam[:,:,40//2])
        exit()
        # dotsInitial = sing.get_singularities(np.angle(beam), axesAll=True, returnDict=False)
        # fig = pl.plot_3D_dots_go(dots, marker={'size': 15, 'color': 'black', 'line': dict(width=185,
        #                                                                                   color='white')})
        # pl.box_set_go(fig, mesh=None, autoDots=dots, perBox=0.05)
        # fig.show()
        # exit()
        dotsDict, dotsInitial = sing.get_singularities(np.angle(beam), axesAll=True, returnDict=True)
        dotsPerfectOnly = [[0, 0, 0]]

        # find_the_group = True
        # if find_the_group

        interpolation_lines = False
        if interpolation_lines:
            dotsInitial = np.array([[41., 57., 92.], [39., 57., 92.4], [43., 57., 91.2], [23., 47., 119.6],
                             [27., 47., 115.2], [25., 45., 122.], [25., 49., 114.], [29., 49., 109.6],
                             [29., 47., 114.4], [27., 49., 111.2], [23., 45., 125.6], [31., 49., 106.8],
                             [25., 47., 117.6], [39., 55., 95.6], [37., 53., 98.4], [35., 55., 96.8],
                             [33., 53., 116.8], [23., 43., 132.8], [25., 41., 145.2], [25., 43., 133.6],
                             [29., 51., 106.4], [31., 53., 121.2], [31., 51., 104.8], [41., 55., 93.6],
                             [33., 51., 103.6], [35., 53., 99.6], [37., 55., 96.4]])

            x_sample = dotsInitial[:, 0]
            y_sample = dotsInitial[:, 1]
            z_sample = dotsInitial[:, 2]
            # import scipy.optimize as optimize
            # params, cov = optimize.curve_fit(fit, np.transpose(XY), Z, guess)
            from scipy import interpolate
            tck = interpolate.splprep([x_sample, y_sample, z_sample], s=2)[0]

            x_knots, y_knots, z_knots = interpolate.splev(tck[0], tck)
            import matplotlib.pyplot as plt
            fig2 = plt.figure(2)
            ax3d = fig2.add_subplot(111, projection='3d')
            ax3d.plot(x_knots, y_knots, z_knots, 'g-')
            ax3d.plot(x_sample, y_sample, z_sample, 'r*')
            fig2.show()
            plt.show()
            exit()

        straight_lines = False
        if straight_lines:
            for dot, value in dotsDict.items():
                neighbours = 0
                xD, yD, zD = dot
                if 1:
                    if dotsDict.get((xD - 1, yD, zD), 10) + dotsDict.get((xD + 1, yD, zD), 10) < 5:
                        dotsPerfectOnly.append(dot)
                    if dotsDict.get((xD, yD - 1, zD), 10) + dotsDict.get((xD, yD + 1, zD), 10) < 5:
                        dotsPerfectOnly.append(dot)
                    if dotsDict.get((xD, yD, zD - 1), 10) + dotsDict.get((xD, yD, zD + 1), 10) < 5:
                        dotsPerfectOnly.append(dot)
                if 1:
                    if dotsDict.get((xD - 1, yD - 1, zD - 1), 10) + dotsDict.get((xD + 1, yD + 1, zD + 1), 10) < 5:
                        dotsPerfectOnly.append(dot)
                    if dotsDict.get((xD + 1, yD - 1, zD - 1), 10) + dotsDict.get((xD - 1, yD + 1, zD + 1), 10) < 5:
                        dotsPerfectOnly.append(dot)
                    if dotsDict.get((xD - 1, yD + 1, zD - 1), 10) + dotsDict.get((xD + 1, yD - 1, zD + 1), 10) < 5:
                        dotsPerfectOnly.append(dot)
                    if dotsDict.get((xD + 1, yD + 1, zD + 1), 10) + dotsDict.get((xD - 1, yD - 1, zD - 1), 10) < 5:
                        dotsPerfectOnly.append(dot)
                if 1:
                    if dotsDict.get((xD - 1, yD, zD - 1), 10) + dotsDict.get((xD + 1, yD, zD + 1), 10) < 5:
                        dotsPerfectOnly.append(dot)
                    if dotsDict.get((xD + 1, yD, zD - 1), 10) + dotsDict.get((xD - 1, yD, zD + 1), 10) < 5:
                        dotsPerfectOnly.append(dot)
                    if dotsDict.get((xD, yD - 1, zD - 1), 10) + dotsDict.get((xD, yD + 1, zD + 1), 10) < 5:
                        dotsPerfectOnly.append(dot)
                    if dotsDict.get((xD, yD + 1, zD - 1), 10) + dotsDict.get((xD, yD - 1, zD + 1), 10) < 5:
                        dotsPerfectOnly.append(dot)
                    if dotsDict.get((xD - 1, yD + 1, zD), 10) + dotsDict.get((xD + 1, yD - 1, zD), 10) < 5:
                        dotsPerfectOnly.append(dot)
                    if dotsDict.get((xD - 1, yD - 1, zD), 10) + dotsDict.get((xD + 1, yD + 1, zD), 10) < 5:
                        dotsPerfectOnly.append(dot)

            # for x in [xD - 1, xD, xD + 1]:
            #     for y in [yD - 1, yD, yD + 1]:
            #         for z in [zD - 1, zD, zD + 1]:
            #             if dotsDict.get((x, y, z), 0) != 0:
            #                 neighbours += 1
            # if neighbours == 2:
            #     dotsPerfectOnly.append(dot)

        dotsGood = np.array(dotsPerfectOnly)

        # print(len(dotsInitial))
        # dotsCutSeparated = sing.dots_filter(dotsInitial, checkValue=np.sqrt(2)*1.01, checkNumber=1)
        # print(len(dotsInitial))
        # dotsCutBigCluster = sing.dots_filter(dotsCutSeparated, checkValue=4, checkNumber=3)
        # dotsFiltered = sing.dots_dens_reduction(dotsCutBigCluster, checkValue=2.5, checkNumber=3)
        # dots = dotsFiltered
        # pl.plot_2D(np.abs(beam[:, :, zRes//2]))
        # pl.plot_scatter_3D(dots[:, 0], dots[:, 1], dots[:, 2], show=True)
        fig = pl.plot_3D_dots_go(dotsInitial, marker={'size': 15, 'color': 'black', 'line': dict(width=185,
                                                                                                 color='white')})
        pl.plot_3D_dots_go(dotsGood, fig=fig, marker={'size': 15, 'color': 'blue', 'line': dict(width=185,
                                                                                                color='white')})
        pl.box_set_go(fig, mesh=None, autoDots=dotsDict, perBox=0.05)
        fig.show()
