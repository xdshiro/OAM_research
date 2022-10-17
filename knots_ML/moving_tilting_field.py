import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np
import scipy.io as sio
import center_beam_search as cbs

if __name__ == '__main__':
    test_hopf_turb_path = 'C:\\WORK\\CODES\\OAM_research\\knots_ML\\' \
                          'data\\Efield_SR_9.000000e-01.mat'
    test_hopf_turb = np.array(sio.loadmat(test_hopf_turb_path, appendmat=False)['Efield'])[:, :, 0]
    # pl.plot_2D(np.abs(test_hopf_turb))
    test_hopf_turb = test_hopf_turb / np.sqrt(np.sum(np.abs(test_hopf_turb)**2))
    shape = np.shape(test_hopf_turb)
    xyMesh = fg.create_mesh_XY(xRes=shape[0], yRes=shape[1])
    # new_xy_mesh = removeShift(xyMesh, 0.5, 0.5)
    beam = cbs.removeTilt(test_hopf_turb, xyMesh, eta=30*np.pi/180, gamma=2*np.pi/180)
    ax = pl.plot_2D(np.abs(beam), x=[-32, 32], y=[-32, 32], show=False)
    pl.plot_scatter_2D(x=3, y= 2.5, ax=ax, color='white', size=150, show=True)
    ax = pl.plot_2D(np.angle(beam), x=[-32, 32], y=[-32, 32], show=False)
    pl.plot_scatter_2D(x=3, y= 2.5, ax=ax, color='black', size=150, show=True)
    exit()

    # ax = pl.plot_2D(np.abs(test_hopf_turb), x=[-32, 32], y=[-32, 32], show=False)
    # pl.plot_scatter_2D(x=2, y=0,ax=ax, color='r', size=100, show=True)
    # exit()
    # plt.scatter()
    pl.plot_2D(np.angle(test_hopf_turb), x=[-32, 32], y=[-32, 32])
    # print(xyMesh)
    # exit()
    # width = cbs.find_width(test_hopf_turb, mesh=xyMesh, width=10, widthStep=1)
    width = 7
    spec = cbs.LG_spectrum(test_hopf_turb,
                           l=(-7, 9), p=(0, 5), mesh=xyMesh, plot=True, width=width)
    # cbs.beam_center_coordinates()
    pl_dict = {'p': (0, 4), 'l': (-5, 3), 'p0': (0, 4), 'l0': (-5, 3)}
    cbs.beamFullCenter(test_hopf_turb, mesh=xyMesh, stepXY=(0.5, 0.5),
                       stepEG=(2 / 180 * np.pi, 1 / 180 * np.pi),
                       **pl_dict, threshold=1,
                       width=width, k0=1, x=0, y=0, eta2=0., gamma=0.)

    plot_trefoil = False
    if plot_trefoil:
        xyzMesh = fg.create_mesh_XYZ(2.1, 2.1, 0.7, 50, 50, 50, zMin=None)
        beam = bp.LG_combination(*xyzMesh,
                                 coefficients=[1.71, -5.66, 6.38, -2.30, -4.36],
                                 modes=[(0, 0), (0, 1), (0, 2), (0, 3), (3, 0)],
                                 width=[1, 1, 1, 1, 1])
        dots = sing.get_singularities(np.angle(beam), axesAll=True)
        # pl.plot_2D(np.abs(beam[:, :, zRes//2]))
        # pl.plot_scatter_3D(dots[:, 0], dots[:, 1], dots[:, 2])
        fig = pl.plot_3D_dots_go(dots)
        pl.box_set_go(fig, mesh=None, autoDots=dots, perBox=0.05)
        fig.show()
