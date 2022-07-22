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
    xMinMax = 3
    yMinMax = 3
    zMinMax = 1
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
    # sing.plot_knot_dots(beam, show=True, axesAll=True)
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
