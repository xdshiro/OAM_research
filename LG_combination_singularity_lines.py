import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np


def plot_line_and_field(beam, xyzMesh):
    beamToPlot = np.abs(beam) / np.abs(beam).max()
    dots = sing.get_singularities(np.angle(beam), bigSingularity=False, axesAll=True)
    dots = fg.dots3D_rescale(dots, xyzMesh)
    fig = pl.plot_3D_density(beamToPlot, mesh=xyzMesh,
                             colorscale='jet',
                             opacityscale=[[0, 0], [0.1, 0], [0.3, 1], [1, 1]], opacity=0.5, surface_count=20,
                             show=False)
    fig.add_trace(pl.plot_3D_dots_go(dots))
    fig.show()


def find_closet_to_point_dot(dots, point):
    ind = 0
    distance = np.inf
    for i, dot in enumerate(dots):
        distanceCurrent = fg.distance_between_points(dot, point)
        if distanceCurrent < distance:
            distance = distanceCurrent
            ind = i
    return ind
if __name__ == '__main__':
    milnor = False
    if milnor:
        xyzMesh = fg.create_mesh_XYZ(5, 5, 2, zMin=None)
        beam = bp.milnor_Pol_u_v_any(xyzMesh, uOrder=4, vOrder=1, H=1)
        dots = sing.get_singularities(np.angle(beam), bigSingularity=False, axesAll=True)
        dots = fg.dots3D_rescale(dots, xyzMesh)
        pl.plot_scatter_3D(dots[:, 0],dots[:, 1],dots[:, 2])
        # plot_line_and_field(beam, xyzMesh)
        exit()
    pyknotid = True
    if pyknotid:
        dotsExp = np.load('C:\\Users\\Dima\\Box\\Knots Exp\\Experimental Data\\dots\\trefoil\\'
                          'Field 3foil noturb\\3foil_noturb_1.npy',  # 25
                          allow_pickle=True).item()
        dots = sing.get_singularities(dotsExp)

        # pl.plot_scatter_3D(dots[:, 0], dots[:, 1], dots[:, 2], size=100)
        dotsKnot = sing.knot_sequence_from_dots(dots, checkValue1=2, checkNumber1=1,
                                           checkValue2=4, checkNumber2=3,
                                           checkValue3=3, checkNumber3=3)
        # pl.plot_scatter_3D(dotsKnot[:, 0], dotsKnot[:, 1], dotsKnot[:, 2], size=100)
        dotsKnot = np.roll(dotsKnot, -1 * find_closet_to_point_dot(dotsKnot, [0, 0, 0]), axis=0)
        sing.plot_knot_pyknotid(dotsKnot, interpolation=400, tube_radius=2.5, tube_points=12, fov=30)
        modesTrefoil = [(0, 0), (0, 1), (0, 2), (0, 3), (3, 0)]
        coeffTrefoil = [1.29, -3.95, 7.49, -3.28, -3.98]
        xMinMax = 2.1
        yMinMax = 2.1
        zMinMax = 0.5
        zRes = 70
        xRes = yRes = 70
        xyzMesh = fg.create_mesh_XYZ(xMinMax, yMinMax, zMinMax, xRes, yRes, zRes, zMin=None)
        beam = bp.LG_combination(*xyzMesh, coefficients=coeffTrefoil, modes=modesTrefoil)
        # sing.plot_knot_dots(beam, show=True, axesAll=False, color='k')
        dots = sing.get_singularities(np.angle(beam), bigSingularity=False, axesAll=True)
        # pl.plot_scatter_3D(dots[:, 0], dots[:, 1], dots[:, 2], size=100)
        dotsKnot = sing.knot_sequence_from_dots(dots, checkValue1=2, checkNumber1=1,
                                                checkValue2=4, checkNumber2=3,
                                                checkValue3=3, checkNumber3=3)[::-1]
        dotsKnot = np.roll(dotsKnot, -1 * find_closet_to_point_dot(dotsKnot, [0, 0, 0]), axis=0)
        # pl.plot_scatter_3D(dotsKnot[:, 0], dotsKnot[:, 1], dotsKnot[:, 2], size=100)
        sing.plot_knot_pyknotid(dotsKnot, interpolation=400, tube_radius=2.5, tube_points=12, clf=True)
        # , antialias=True,
        #                                 light_dir=(180,90,-50)
        # dots = fg.dots3D_rescale(dots, xyzMesh)
        exit()
    xMinMax = 2.1
    yMinMax = 2.1
    zMinMax = 0.5
    zRes = 70
    xRes = yRes = 70
    xyzMesh = fg.create_mesh_XYZ(xMinMax, yMinMax, zMinMax, xRes, yRes, zRes, zMin=None)
    # modes = [(0, 0), (0, 1), (0, 2), (1, 0), (2, 0), (2, 1), (2, 2)]
    # modes = [(0, 0), (0, 1), (1, 0), (1, 1)]
    # lpArray = [0, 0, 0, 0]
    # diapason = [1, 1, 1, 1]
    modes = [(0, 0), (0, 1), (1, 0), (1, 1), (-1, 0), (-1, 1)]
    lpArray = [0, 0, 0, 0, 0, 0]
    lpDiapason = [1, 1, 1, 1, 1, 1]
    phaseArray = [0, 0, 0, 0, 0, 0]
    phaseDiapason = [0, 0, 0, 0, 0, 0]
    wArray = [1, 1, 1, 1, 1, 1]
    wDiapason = [0, 0, 0, 0, 0, 0]
    modesHopf = [(0, 0), (0, 1), (0, 2), (2, 0)]
    coeffHopf = [2.96, -6.23, 4.75, -5.49]
    modesTrefoil = [(0, 0), (0, 1), (0, 2), (0, 3), (3, 0)]
    coeffTrefoil = [1.29, -3.95, 7.49, -3.28, -3.98]
    beam = bp.LG_combination(*xyzMesh, coefficients=coeffHopf, modes=modesHopf)
    # sing.plot_knot_dots(beam, show=True, axesAll=False, color='k')
    dots = sing.get_singularities(np.angle(beam), bigSingularity=False, axesAll=True)
    dots = fg.dots3D_rescale(dots, xyzMesh)
    # dots = sing.get_singularities(dots)
    fig = pl.plot_3D_dots_go(dots)
    # xyzMesh = fg.create_mesh_XYZ(xMinMax, yMinMax, zMinMax, xRes, yRes, zRes, zMin=None)
    # dotsExp = np.load(f'C:\\Users\\Dima\\Box\\Knots Exp\\Experimental Data\\dots'
    #                   f'\\trefoil\\Field 3foil noturb\\3foil_noturb_1.npy',
    #                allow_pickle=True).item()
     # dots = sing.get_singularities(dotsExp)
    # xMinMax = 200
    # yMinMax = 200
    # zMinMax = 80
    # fig.update_layout(
    #     scene=dict(
    #         aspectratio_x=1.95, aspectratio_y=1.95, aspectratio_z=1,
    #         xaxis=dict(range=[70, 150]),
    #         yaxis=dict(range=[70, 150]),
    #         zaxis=dict(range=[10, 55]), ),
    #     # width=700,
    #     # margin=dict(r=0, l=10, b=70, t=10)
    # )
    # print(dots)
    # exit()

    # fig.update_layout(title='Mt Bruno Elevation', autosize=True)
    # fig.update_layout(yaxis_range=[-4,4])
    # import plotly.io as pio
    # pio.renderers.default = "png"
    # fig.layout.xaxis.range = (0,1)
    # iplot(fig)
    fig.update_layout(
        scene=dict(
            aspectratio_x=2, aspectratio_y=2, aspectratio_z=1,
            xaxis=dict(range=[-xMinMax, xMinMax]),
            yaxis=dict(range=[-yMinMax, yMinMax]),
            zaxis=dict(range=[-zMinMax, zMinMax]),),
        # width=700,
        # margin=dict(r=0, l=10, b=70, t=10)
    )
    fig.show()
    exit()

    while False:
        coeff = fg.random_list(lpArray, lpDiapason, diapason_complex=None)
        phase = fg.random_list(phaseArray, phaseDiapason, diapason_complex=None)
        width = fg.random_list(wArray, wDiapason, diapason_complex=None)

        coeffComp = [c * np.exp(1j * p) for c, p in zip(coeff, phase)]
        # beam = bp.LG_combination(*xyzMesh, coefficients=coeff, modes=modes)
        # pl.plot_2D(np.abs(beam[:, :, zRes // 2]))
        # sing.plot_knot_dots(beam, show=True, axesAll=True)
        beam = bp.LG_combination(*xyzMesh, coefficients=coeffComp, modes=modes, width=width)
        # sing.plot_knot_dots(beam, show=True, axesAll=True)
        beamToPlot = np.abs(beam) / np.abs(beam).max()
        dots = sing.get_singularities(np.angle(beam), bigSingularity=False, axesAll=True)
        fig = pl.plot_3D_density(beamToPlot, colorscale='jet',
                                 opacityscale=[[0, 0], [0.1, 0], [0.3, 1], [1, 1]], opacity=0.5, surface_count=20,
                                 show=False)
        fig.add_trace(pl.plot_3D_dots_go(dots))
        coeff = [round(c, 3) for c in coeff]
        phase = [round(p, 3) for p in phase]
        width = [round(w, 3) for w in width]
        print(f'{modes}: amp: {coeff}, phase:{phase}2pi, width: {width}')
        fig.layout.title.text = f'{modes}: amp: {coeff}, phase:{phase}2pi, width: {width}'
        fig.show()

    coeff = [1, 1, 2, 1, 2, 1]
    phase = [0.1, 0, 0.3, 0, -0.1, 0]
    width = [2, 1, 1.2, 1, 2, 1]
    coeffComp = [c * np.exp(1j * p * 2 * np.pi) for c, p in zip(coeff, phase)]
    beam = bp.LG_combination(*xyzMesh, coefficients=coeffComp, modes=modes, width=width)
    beamToPlot = np.abs(beam) / np.abs(beam).max()
    dots = sing.get_singularities(np.angle(beam), bigSingularity=False, axesAll=True)
    xyzMesh = fg.create_mesh_XYZ(xMinMax, yMinMax, zMinMax * 3, xRes, yRes, zRes, zMin=None)
    dots = fg.dots3D_rescale(dots, xyzMesh)
    fig = pl.plot_3D_density(beamToPlot, mesh=xyzMesh,
                             colorscale='jet',
                             opacityscale=[[0, 0], [0.1, 0], [0.3, 1], [1, 1]], opacity=0.4, surface_count=20,
                             show=False)
    fig.add_trace(pl.plot_3D_dots_go(dots))
    coeff = [round(c, 3) for c in coeff]
    phase = [round(p, 3) for p in phase]
    width = [round(w, 3) for w in width]
    fig.layout.title.text = f'{modes}: amp: {coeff}, phase:{phase}2pi, width: {width}'
    fig.show()

exit()
# pl.plot_3D_density(beamToPlot, colorscale='jet',
#                    opacity=0.0, surface_count=20, show=True)
