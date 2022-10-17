import my_functions.functions_general as fg
import my_functions.singularities as sing
import my_functions.plotings as pl
import my_functions.beams_and_pulses as bp
import numpy as np
from PIL import Image
import os


def get_dots_from_pic(fileName, size=(100, 100), res=40, valueRed=None, plot=False, cutCenter=False):
    from PIL import Image
    image = Image.open(fileName)
    image = image.resize(size, Image.ANTIALIAS)
    # convert image to numpy array
    data = np.asarray(image, dtype=int)
    if cutCenter:
        print(np.shape(data))
        dataInt = np.sqrt(data[:, :, 0] ** 2 + data[:, :, 1] ** 2 + data[:, :, 2] ** 2)
        centerX, centerY = find_center_mass_2D(dataInt, integer=True)
        data = data[centerX - res: centerX + res + 1, centerY - res: centerY + res + 1, :]
    # pl.plot_2D(data[:, :, 1])
    # exit()
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


def find_center_mass_2D(array, integer=False):
    shape = np.shape(array)
    array = np.copy(array / np.abs(array).max())
    sumMrX = 0
    sumMrY = 0
    sumM = 0
    for i in range(shape[0]):
        for j in range(shape[1]):
            sumMrX += i * array[i, j]
            sumMrY += j * array[i, j]
            sumM += array[i, j]
    centerX = sumMrX / sumM
    centerY = sumMrY / sumM
    if integer:
        centerX, centerY = int(centerX), int(centerY)
    return centerX, centerY


def dots_from_directory(directoryName, **kwargs):
    listOfFiles = [f for f in os.listdir(directoryName) if f.endswith('png')]
    dots3D = []
    for file in listOfFiles:
        print(file)
        dots = get_dots_from_pic(directoryName + '\\' + file, **kwargs)
        z = float(file[:-4])
        dotsWithZ = add_Z_to_dots(dots, z)
        for dot in dotsWithZ:
            dots3D.append(dot)
    return np.array(dots3D)


def cut_all_plots_from_directory(directoryName, res):
    from PIL import Image
    listOfFiles = [f for f in os.listdir(directoryName) if f.endswith('.png')]
    print(listOfFiles)
    for file in listOfFiles:
        image = Image.open(directoryName + file)
        # image = image.resize(size, Image.ANTIALIAS)
        # convert image to numpy array
        data = np.asarray(image, dtype=int)
        print(np.shape(data))
        # dataInt = np.sqrt(data[:, :, 0] ** 2 + data[:, :, 1] ** 2 + data[:, :, 2] ** 2)
        centerX, centerY = find_center_mass_2D(data, integer=True)
        print(file, centerX, centerY)
        # pl.plot_2D(data)
        # exit()
        data = data[centerX - res: centerX + res + 1, centerY - res: centerY + res + 1]
        Image.fromarray(data).save(file[:-4] + 'res.png')

def dots_recentering(dots):
    vals, indxStart, count = np.unique(dots[:, 2], return_counts=True, return_index=True)
    # print(dots)
    xAvgTot = 0
    yAvgTot = 0
    countPlanes = 0
    for i, value in enumerate(count):
        if value == 6:
            values = dots[indxStart[i]:indxStart[i] + 6]
            xAvgTot += sum(values[:, 0]) / 6
            yAvgTot += sum(values[:, 1]) / 6
            countPlanes += 1
    xTot = xAvgTot / countPlanes
    yTot = yAvgTot / countPlanes
    for i, value in enumerate(count):
        if value == 6:
            values = dots[indxStart[i]:indxStart[i] + 6]
            xAvg = sum(values[:, 0]) / 6
            yAvg = sum(values[:, 1]) / 6
            # print(dots[indxStart[i]:indxStart[i] + 6][:, 0])
            dots[indxStart[i]:indxStart[i] + 6][:, 0] += xTot - xAvg
            dots[indxStart[i]:indxStart[i] + 6][:, 1] += yTot - yAvg
    return dots
if __name__ == '__main__':
    from PIL import Image

    # cut_all_plots_from_directory('.\\trefoil_9_2\\', res=42)
    # exit()
    # image = Image.open('.\\DATA\\Jiannan\\1.jpg')
    # # image = image.resize(size, Image.ANTIALIAS)
    # # convert image to numpy array
    # data = np.asarray(image, dtype=int)
    # dataInt = np.sqrt(data[:, :, 0] ** 2 + data[:, :, 1] ** 2 + data[:, :, 2] ** 2)
    #
    # centerX, centerY = find_center_mass_2D(dataInt, integer=True)
    # res = 150
    # data = data[centerX - res: centerX + res + 1, centerY - res: centerY + res + 1, :]
    # # ax= pl.plot_2D(dataInt, show=True)
    # ax = pl.plot_2D(data[:, :, 1], show=True)
    # # pl.plot_scatter_2D(centerX, centerY, size=300, ax=ax, color='r')
    #
    # print(centerX, centerY)
    # exit()
    # dots = get_dots_from_pic('.\\DATA\\Jiannan\\1.jpg')
    # pl.plot_scatter_2D(dots[:, 0], dots[:, 1],
    #                    size=100, color='r', xlim=[0, 100], ylim=[0, 100], axis_equal=True)
    # exit()
    dots = dots_from_directory('.\\10_10_2022\\dots_cropped', plot=False, cutCenter=False)
    for i in range(len(dots)):
        if dots[i][2] < 0.5:

            dots[i][0] -= 40

            dots[i][0] *= 0.8

            dots[i][0] += 40

            dots[i][1] -= 40
            dots[i][1] *= 0.8
            dots[i][1] += 40
            # dots[i][1], dots[i][0] = dots[i][0], dots[i][1]

    # exit()
    # dots = dots_recentering(dots)
    # print(dots)
    #
    # print(dots[:, 2])

    # print(dots)

    fig = pl.plot_3D_dots_go(dots, show=False, marker={'size': 20, 'color': 'black', 'line': dict(width=185,
                                                                                                  color='white')})
    pl.box_set_go(fig, mesh=None, autoDots=dots, perBox=0.05)
    fig.show()
    exit()
