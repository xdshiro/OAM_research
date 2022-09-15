import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np

if __name__ == '__main__':
    'C:\\WORK\\CODES\\OAM_research\\Wenhao_hypergrating\data'
    with open('C:\\WORK\\CODES\\OAM_research\\Wenhao_hypergrating\\data\\Field_xy_1.npy', 'rb') as f:
        a = np.load(f)
        b = np.load(f)

    pl.plot_2D(np.abs(a), show=False)
    pl.plot_2D(np.angle(a), show=False)
    pl.plot_2D(np.abs(b), show=False)
    pl.plot_2D(np.angle(b))

    Jza = sing.Jz_calc_no_conj(a)
    Jzb = sing.Jz_calc_no_conj(b)
    Wa = sing.W_energy(a)
    Wb = sing.W_energy(b)

    print('total combined: ',(Jza + Jzb) / (Wa + Wb))
    # print(np.abs(b).max())