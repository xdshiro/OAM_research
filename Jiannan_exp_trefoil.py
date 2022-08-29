import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np
from PIL import Image
import os


def get_dots_from_pic(fileName, size=(100, 100), valueRed=None, plot=False):
    from PIL import Image
    image = Image.open(fileName)
    image = image.resize(size, Image.ANTIALIAS)
    # convert image to numpy array
    data = np.asarray(image, dtype=int)
    dataRed = data[:, :, 0] - data[:, :, 1] - data[:, :, 2]
    # pl.plot_2D(data[:, :, 0])
    # pl.plot_2D(data[:, :, 1])
    # pl.plot_2D(data[:, :, 2])
    # pl.plot_2D(dataRedOnly[:, :])
    if valueRed is None:
        max = dataRed.max()
        valueRed = 0.6 * max
    dataRed[dataRed < valueRed] = 0
    if plot:
        pl.plot_2D(dataRed)
    shape = np.shape(dataRed)
    dots = []
    circle = 2
    for i in range(shape[0]):
        for j in range(shape[1]):
            if dataRed[i, j] > 0:
                dataRed[i - circle:i + circle + 1, j - circle:j + circle + 1] = 0
                dataRed[i, j] = 1
                dots.append([i, j])
    dots = np.array(dots)
    return dots


def add_Z_to_dots(dots, z):
    zArray = np.array([np.zeros(len(dots)) + z])
    dotsWithZ = np.concatenate((dots, zArray.T), axis=1)
    return dotsWithZ


def dots_from_directory(directoryName):
    listOfFiles = [f for f in os.listdir(directoryName) if f.endswith('jpg')]
    dots3D = []
    for file in listOfFiles:
        print(file)
        dots = get_dots_from_pic(directoryName + '\\' + file)
        z = float(file[:-4])
        dotsWithZ = add_Z_to_dots(dots, z)
        for dot in dotsWithZ:
            dots3D.append(dot)
    return np.array(dots3D)


if __name__ == '__main__':
    from PIL import Image

    image = Image.open('.\\DATA\\Jiannan\\1.jpg')
    # image = image.resize(size, Image.ANTIALIAS)
    # convert image to numpy array
    data = np.asarray(image, dtype=int)
    dataInt = np.sqrt(data[:, :, 0] ** 2 + data[:, :, 1] ** 2 + data[:, :, 2] ** 2)
    # pl.plot_2D(dataInt)
    shape = np.shape(dataInt)
    sumMr = 0
    sumM = 0
    for i in range(shape[0]):
        for j in range(shape[1]):

    exit()
    dots = get_dots_from_pic('.\\DATA\\Jiannan\\1.jpg')
    pl.plot_scatter_2D(dots[:, 0], dots[:, 1],
                       size=100, color='r', xlim=[0, 100], ylim=[0, 100], axis_equal=True)
    exit()
    dots = dots_from_directory('.\\DATA\\Jiannan')
    pl.plot_3D_dots_go(dots, show=True)
    exit()
