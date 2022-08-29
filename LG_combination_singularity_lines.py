import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np

def plot_line_and_field(beam, xyzMesh, show=False, opacityscale=([0, 0], [0.2, 0], [0.3, 0.5], [1, 1]), opacity=0.4):
    beamToPlot = np.abs(beam) / np.abs(beam).max()
    dots = sing.get_singularities(np.angle(beam), bigSingularity=False, axesAll=True)
    dots = fg.dots3D_rescale(dots, xyzMesh)
    fig = pl.plot_3D_density(beamToPlot, mesh=xyzMesh,
                             colorscale='jet',
                             opacityscale=opacityscale, opacity=opacity, surface_count=20,
                             show=False)
    pl.plot_3D_dots_go(dots, marker={'size': 8, 'color': 'black'}, fig=fig)
    if show:
        fig.show()
    return fig



def find_closet_to_point_dot(dots, point):
    ind = 0
    distance = np.inf
    for i, dot in enumerate(dots):
        distanceCurrent = fg.distance_between_points(dot, point)
        if distanceCurrent < distance:
            distance = distanceCurrent
            ind = i
    return ind




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
        #
        # dotsExp = np.load('C:\\Users\\Dima\\Box\\Knots Exp\\Experimental Data\\dots\\trefoil\\'
        #                   'Field 3foil noturb\\3foil_noturb_1.npy',  # 25
        #                   allow_pickle=True).item()
        # dots = sing.get_singularities(dotsExp)
        #
        #
        # # pl.plot_scatter_3D(dots[:, 0], dots[:, 1], dots[:, 2], size=100)
        # dotsKnot = sing.knot_sequence_from_dots(dots, checkValue1=2, checkNumber1=1,
        #                                         checkValue2=4, checkNumber2=3,
        #                                         checkValue3=3, checkNumber3=3)
        # # pl.plot_scatter_3D(dotsKnot[:, 0], dotsKnot[:, 1], dotsKnot[:, 2], size=100)
        # dotsKnot = np.roll(dotsKnot, -1 * find_closet_to_point_dot(dotsKnot, [0, 0, 0]), axis=0)
        #
        #
        #
        #
        # # dotsKnot = dots_move_center(dotsKnot)
        # # fig = pl.plot_3D_dots_go(dotsKnot)
        # # fig.show()
        # # print(dotsKnot)
        # # exit()
        # sing.plot_knot_pyknotid(dotsKnot, interpolation=300, tube_radius=2.5, per=True, add_closure=False,
        #                         tube_points=14, fov=0, flip=(False, False, True))
        # plot_vispy_tube(dotsKnot[20:])

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
        sing.plot_knot_pyknotid(dotsKnot, interpolation=300, tube_radius=2.0, per=False, add_closure=False,
                                tube_points=14, fov=0, flip=(False, False, True))
        # , antialias=True,
        #                                 light_dir=(180,90,-50)
        # dots = fg.dots3D_rescale(dots, xyzMesh)
        exit()

    random_combination = False
    while random_combination:
        xMinMax = 2.1
        yMinMax = 2.1
        zMinMax = 0.5
        zRes = 80
        xRes = yRes = 80
        xyzMesh = fg.create_mesh_XYZ(xMinMax, yMinMax, zMinMax, xRes, yRes, zRes, zMin=None)
        modes = [(0, 0), (0, 1), (1, 0), (1, 1), (-1, 0), (-1, 1)]
        lpArray = [0.5] * 5
        lpDiapason = [0.5] * 5
        phaseArray = [np.pi] * 6
        phaseDiapason = [np.pi] * 6
        wArray = [0.5] * 5
        wDiapason = [0.5] * 5
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
        # fig = pl.plot_3D_density(beamToPlot, colorscale='jet',
        #                          opacityscale=[[0, 0], [0.1, 0], [0.3, 1], [1, 1]], opacity=0.5, surface_count=20,
        #                          show=False)
        fig = pl.plot_3D_dots_go(dots)
        coeff = [round(c, 3) for c in coeff]
        phase = [round(p, 3) for p in phase]
        width = [round(w, 3) for w in width]
        print(f'{modes}: amp: {coeff}, phase:{phase}2pi, width: {width}')
        fig.layout.title.text = f'{modes}: amp: {coeff}, phase:{phase}2pi, width: {width}'
        fig.show()
    single_plot = False
    if single_plot:
        magnification = 1
        xMinMax = 2.1 / magnification
        yMinMax = 2.1 / magnification
        zMinMax = 0.5 / magnification
        zRes = 50
        xRes = yRes = 50
        xyzMesh = fg.create_mesh_XYZ(xMinMax, yMinMax, zMinMax, xRes, yRes, zRes, zMin=None)
        modes = [(0, 0), (0, 1), (0, 2), (1, 0), (2, 0), (2, 1), (2, 2)]
        coeff = [1, 1, 2, 1, 2, 1]
        phase = [0.1, 0, 0.3, 0, -0.1, 0]
        width = [2, 1, 1.2, 1, 2, 1]
        coeffComp = [c * np.exp(1j * p * 2 * np.pi) for c, p in zip(coeff, phase)]
        # -1 + w^2, -w^2, -2 w
        modes = [(0, 0), (0, 1), (3, 0), (1,0)]
        w = 2
        coeffComp = [1 + w ** 2, -w ** 2, -2 * w * 2]
        # coeffComp = [1, 1, 0, 0]
        # width = [1, 1, 0.5, 1]
        beam = bp.LG_combination(*xyzMesh, coefficients=coeffComp, modes=modes, width=width)
        # beam *= np.exp(1j * 2 * fg.phi(xyzMesh[0], xyzMesh[1]))
        pl.plot_2D(abs(beam[:,:,zRes//2]))
        # pl.plot_2D(np.angle(beam[:,:,zRes//2]))
        pl.plot_2D(abs(beam[xRes//2,:,:]))
        # pl.plot_2D(np.angle(beam[xRes//2,:,:]))
        # exit()
        # pl.plot_2D(abs(beam[:,:,zRes//2 + 20]))
        # pl.plot_2D(abs(beam[:,:,zRes//2 - 30]))
        # exit()
        # exit()
        # dots = sing.get_singularities(np.angle(beam), bigSingularity=False, axesAll=True)
        # dots = fg.dots3D_rescale(dots, xyzMesh)
        fig = plot_line_and_field(beam, xyzMesh, opacityscale=([0, 1], [0.2, 1], [0.3, 0], [1, 0]), opacity=0.3)
        coeff = [round(c, 3) for c in coeff]
        phase = [round(p, 3) for p in phase]
        width = [round(w, 3) for w in width]
        fig.update_layout(
            scene=dict(
                aspectratio_x=2, aspectratio_y=2, aspectratio_z=2,
                xaxis=dict(range=[-xMinMax, xMinMax]),
                yaxis=dict(range=[-yMinMax, yMinMax]),
                zaxis=dict(range=[-zMinMax, zMinMax]), ),
            # width=700,
            # margin=dict(r=0, l=10, b=70, t=10)
        )
        fig.layout.title.text = f'{modes}: amp: {coeff}, phase:{phase}2pi, width: {width}'
        fig.show()

    plot_normal_hopf = True
    if plot_normal_hopf:

        magnification = 1
        xMinMax = 2.1 / magnification
        yMinMax = 2.1 / magnification
        zMinMax = 0.7 / magnification
        zRes = 50
        xRes = yRes = 50
        xyzMesh = fg.create_mesh_XYZ(xMinMax, yMinMax, zMinMax, xRes, yRes, zRes, zMin=None)
        modes = [(0, 0), (0, 1), (0, 2), (2, 0)]
        coeff = [2.63, -6.32, 4.21, -5.95]
        modes = [(0, 0), (0, 1), (0, 2), (0, 3), (3, 0)]
        coeff = [1.71, -5.66, 6.38, -2.30, -4.36]
        width = [1, 1, 1, 1, 1]
        beam = bp.LG_combination(*xyzMesh, coefficients=coeff, modes=modes, width=width)
        # pl.plot_2D(np.abs(beam[:, :, zRes//2]))
        plot_line_and_field(beam, xyzMesh, show=True)
