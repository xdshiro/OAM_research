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


def dots_move_center(dots):
    """
    moving dots to the center of the object
    """
    center = np.sum(dots, axis=0) / len(dots)
    return dots - center

def plot_vispy_tube(dots, clf=False, tube_radius=40, shading='smooth', fov=20):
        from vispy import scene
        from colorsys import hsv_to_rgb
        # def clear_vispy_canvas():
        #     global vispy_canvas
        #     if vispy_canvas is None:
        #         return
        #     try:
        #         vispy_canvas.unfreeze()
        #     except AttributeError:  # depends on Vispy version
        #         pass
        #     vispy_canvas.central_widget.remove_widget(vispy_canvas.view)
        #     vispy_canvas.view = vispy_canvas.central_widget.add_view()
        #
        # if clf:
        #     vispy_canvas.central_widget.remove_widget(vispy_canvas.view)
        #     vispy_canvas.view = vispy_canvas.central_widget.add_view()
        canvas = scene.SceneCanvas(keys='interactive')
        view = canvas.central_widget.add_view()
        colors = np.linspace(0, 1, len(dots))
        colors = np.array([hsv_to_rgb(c, 1, 1) for c in colors])
        l = scene.visuals.Tube(dots,
                               color=colors,
                               shading=shading,
                               tube_points=8,
                               radius=tube_radius,
                               closed=True)
        # if zero_centroid:
        #     l.transform = MatrixTransform()
        #     # l.transform = scene.transforms.AffineTransform()
        #     l.transform.translate(-1 * n.average(points, axis=0))
        # view.add(l)

        flip = (True, True, False)
        view.camera = scene.TurntableCamera(fov=fov, flip=flip, distance=7.5 * np.max(np.abs(dots)))
        view.camera.set_range((-20, 20), (-20, 20), (-20, 20))
        canvas.show()
        canvas.app.run()


# def plot_line_vispy(points, clf=True, tube_radius=1.,
#                     colour=None, zero_centroid=True,
#                     closed=True, mus=None,
#                     cmap=None, fov=0, flip=(False, False, False),
#                     tube_points=8, **kwargs):
#     # Add an extra point to fix tube drawing bug
#     last_tangent = points[-1] - points[-2]
#     points = n.vstack([points, points[-1] + 0.0001 * last_tangent])
#
#     ensure_vispy_canvas()
#     if clf:
#         clear_vispy_canvas()
#     canvas = vispy_canvas
#     from vispy import app, scene, color
#
#     if isinstance(cmap, str):
#         from matplotlib.cm import get_cmap
#         mpl_cmap = get_cmap(cmap)
#         cmap = lambda v: n.array(mpl_cmap(v))
#     cmap = cmap or (lambda c: hsv_to_rgb(c, 1, 1))
#
#     if colour is None:
#         colours = n.linspace(0, 1, len(points))
#         colours = n.array([cmap(c) for c in colours])
#     else:
#         colours = color.ColorArray(colour)
#
#     if mus is not None:
#         colours = n.array([hsv_to_rgb(c, 1, 1) for c in mus])
#
#     l = scene.visuals.Tube(points, color=colours,
#                            shading='smooth',
#                            radius=tube_radius,
#                            tube_points=tube_points,
#                            closed=closed)
#
#     canvas.view.add(l)
#     # canvas.view.camera = 'arcball'
#     canvas.view.camera = scene.ArcballCamera(fov=30, flip=flip, distance=7.5 * n.max(
#         n.abs(points)))
#     # canvas.view.camera = scene.TurntableCamera(fov=30)
#     if zero_centroid:
#         l.transform = MatrixTransform()
#         # l.transform = scene.transforms.AffineTransform()
#         l.transform.translate(-1 * n.average(points, axis=0))
#
#     canvas.show()
#     # import ipdb
#     # ipdb.set_trace()
#     return canvas

if __name__ == '__main__':
    milnor = False
    if milnor:
        xyzMesh = fg.create_mesh_XYZ(5, 5, 2, zMin=None)
        beam = bp.milnor_Pol_u_v_any(xyzMesh, uOrder=4, vOrder=1, H=1)
        dots = sing.get_singularities(np.angle(beam), bigSingularity=False, axesAll=True)
        dots = fg.dots3D_rescale(dots, xyzMesh)
        pl.plot_scatter_3D(dots[:, 0], dots[:, 1], dots[:, 2])
        # plot_line_and_field(beam, xyzMesh)
        exit()
    pyknotid = False
    if pyknotid:
        # def plot_vispy_tube(points, clf=True, tube_radius=1.,
        #                     colour=None, #zero_centroid=True,
        #                     closed=False, mus=None,
        #                     cmap=None,
        #                     tube_points=8, **kwargs):
        #     from colorsys import hsv_to_rgb
        #     def clear_vispy_canvas():
        #         global vispy_canvas
        #         if vispy_canvas is None:
        #             return
        #         try:
        #             vispy_canvas.unfreeze()
        #         except AttributeError:  # depends on Vispy version
        #             pass
        #         vispy_canvas.central_widget.remove_widget(vispy_canvas.view)
        #         vispy_canvas.view = vispy_canvas.central_widget.add_view()
        #
        #     def ensure_vispy_canvas():
        #         global vispy_canvas
        #         if vispy_canvas is None:
        #             from vispy import app, scene
        #             canvas = scene.SceneCanvas(keys='interactive', bgcolor='white')
        #             try:
        #                 canvas.unfreeze()
        #             except AttributeError:  # depends on Vispy version
        #                 pass
        #             canvas.view = canvas.central_widget.add_view()
        #             vispy_canvas = canvas
        #     # Add an extra point to fix tube drawing bug
        #     # last_tangent = points[-1] - points[-2]
        #     # points = n.vstack([points, points[-1] + 0.0001 * last_tangent])
        #
        #     ensure_vispy_canvas()
        #     if clf:
        #         clear_vispy_canvas()
        #     canvas = vispy_canvas
        #     from vispy import app, scene, color
        #
        #     if isinstance(cmap, str):
        #         from matplotlib.cm import get_cmap
        #         mpl_cmap = get_cmap(cmap)
        #         cmap = lambda v: np.array(mpl_cmap(v))
        #     cmap = cmap or (lambda c: hsv_to_rgb(c, 1, 1))
        #
        #     if colour is None:
        #         colours = np.linspace(0, 1, len(points))
        #         colours = np.array([cmap(c) for c in colours])
        #     else:
        #         colours = color.ColorArray(colour)
        #
        #     if mus is not None:
        #         colours = np.array([hsv_to_rgb(c, 1, 1) for c in mus])
        #
        #     l = scene.visuals.Tube(points, color=colours,
        #                            shading='smooth',
        #                            radius=tube_radius,
        #                            tube_points=tube_points,
        #                            closed=closed)
        #
        #     canvas.view.add(l)
        #     # canvas.view.camera = 'arcball'
        #     canvas.view.camera = scene.ArcballCamera(fov=30, distance=7.5 * np.max(
        #         np.abs(points)))
        #     # canvas.view.camera = scene.TurntableCamera(fov=30)
        #     # if zero_centroid:
        #     #     l.transform = MatrixTransform()
        #     #     # l.transform = scene.transforms.AffineTransform()
        #     #     l.transform.translate(-1 * n.average(points, axis=0))
        #
        #     canvas.show()
        #     # import ipdb
        #     # ipdb.set_trace()
        #     return canvas

        # dotsExp = np.load('C:\\Users\\Dima\\Box\\Knots Exp\\Experimental Data\\dots\\trefoil\\'
        #                   'Field 3foil noturb\\3foil_noturb_1.npy',  # 25
        #                   allow_pickle=True).item()

        dotsExp = np.load('C:\\Users\\Cmex-\\Box\\Knots Exp\\Experimental Data\\dots\\trefoil\\'
                          'Field 3foil noturb\\3foil_noturb_1.npy',  # 25
                          allow_pickle=True).item()
        dots = sing.get_singularities(dotsExp)


        # pl.plot_scatter_3D(dots[:, 0], dots[:, 1], dots[:, 2], size=100)
        dotsKnot = sing.knot_sequence_from_dots(dots, checkValue1=2, checkNumber1=1,
                                                checkValue2=4, checkNumber2=3,
                                                checkValue3=3, checkNumber3=3)
        # pl.plot_scatter_3D(dotsKnot[:, 0], dotsKnot[:, 1], dotsKnot[:, 2], size=100)
        dotsKnot = np.roll(dotsKnot, -1 * find_closet_to_point_dot(dotsKnot, [0, 0, 0]), axis=0)




        # dotsKnot = dots_move_center(dotsKnot)
        # fig = pl.plot_3D_dots_go(dotsKnot)
        # fig.show()
        # print(dotsKnot)
        # exit()
        sing.plot_knot_pyknotid(dotsKnot, interpolation=300, tube_radius=2.5, per=True, add_closure=False,
                                tube_points=14, fov=0, flip=(False, False, True))
        plot_vispy_tube(dotsKnot)
        exit()
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
        dotsKnot = np.roll(dotsKnot, -1 * find_closet_to_point_dot(dotsKnot, [0, 0, 0]) + 40, axis=0)
        # pl.plot_scatter_3D(dotsKnot[:, 0], dotsKnot[:, 1], dotsKnot[:, 2], size=100)
        sing.plot_knot_pyknotid(dotsKnot, interpolation=300, tube_radius=2.5, per=False, add_closure=False,
                                tube_points=14, fov=0, flip=(False, False, True))
        # , antialias=True,
        #                                 light_dir=(180,90,-50)
        # dots = fg.dots3D_rescale(dots, xyzMesh)
        exit()
    xMinMax = 2.1
    yMinMax = 2.1
    zMinMax = 0.5
    zRes = 200
    xRes = yRes = 200
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
    beam = bp.LG_combination(*xyzMesh, coefficients=coeffTrefoil, modes=modesTrefoil)
    xyz = fg.arrays_from_mesh(xyzMesh)
    # pl.plot_2D(np.angle(beam[:,:,int(170 * 1.5)]), x=xyz[0], y=xyz[1], map='jet')
    #
    # exit()
    # sing.plot_knot_dots(beam, show=True, axesAll=False, color='k')
    dots = sing.get_singularities(np.angle(beam), bigSingularity=False, axesAll=True)
    dots = fg.dots3D_rescale(dots, xyzMesh)
    # dots = sing.get_singularities(dots)
    fig = pl.plot_3D_dots_go(dots, marker={'size':15, 'color':'black'})
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
            aspectratio_x=2, aspectratio_y=2, aspectratio_z=2,
            xaxis=dict(range=[-xMinMax, xMinMax]),
            yaxis=dict(range=[-yMinMax, yMinMax]),
            zaxis=dict(range=[-zMinMax, zMinMax]), ),
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
