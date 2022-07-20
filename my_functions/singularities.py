"""
This module has classes of different singularities and functions processing singularities
"""
# import functions_OAM_knots as fOAM
# import functions_general as fg
# import matplotlib.pyplot as plt
import numpy as np
# import functions_general as fg
import pyknotid.spacecurves as sp
import my_functions.beams_and_pulses as bp
import my_functions.plotings as pl
import my_functions.functions_general as fg
import matplotlib.pyplot as plt
from scipy import integrate

"""
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
We can make a graph for tsp, so it is not searching for all the dots, only close z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
trefoilW16.fill_dotsList() this function is what makes everything slow 

Stop if the distance is too big. To find a hopf 
"""


def plot_knot_dots(field, bigSingularity=False, axesAll=True,
                   size=plt.rcParams['lines.markersize'] ** 2, color=None, show=False):
    """
    ploting the 3d scatters (or 2d) from the field or from the dict with dots
    :param field: can be complex field or dictionary with dots to plot
    :param bigSingularity: cut_non_oam
    :param axesAll: cut_non_oam
    :param size: dots size
    :param color: dots color
    :return:
    """
    if isinstance(field, dict):
        dotsOnly = field
    else:
        dotsFull, dotsOnly = cut_non_oam(np.angle(field),
                                         bigSingularity=bigSingularity, axesAll=axesAll)
    dotsPlus = np.array([list(dots) for (dots, OAM) in dotsOnly.items() if OAM == 1])
    dotsMinus = np.array([list(dots) for (dots, OAM) in dotsOnly.items() if OAM == -1])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if len(np.shape(dotsPlus)) == 2:
        pl.plot_scatter_3D(dotsPlus[:, 0], dotsPlus[:, 1], dotsPlus[:, 2], ax=ax, size=size, color=color)
        if len(np.shape(dotsMinus)) == 2:
            pl.plot_scatter_3D(dotsMinus[:, 0], dotsMinus[:, 1], dotsMinus[:, 2], ax=ax, size=size, color=color)
    else:
        if len(np.shape(dotsPlus)) == 2:
            pl.plot_scatter_3D(dotsMinus[:, 0], dotsMinus[:, 1], dotsMinus[:, 2], ax=ax, size=size, color=color)
        else:
            print(f'no singularities to plot')
    if show:
        plt.show()
    return ax


def plane_singularities_finder_9dots(E, circle, value, nonValue, bigSingularity):
    """
    cut_non_oam helper. see cut_non_oam for more details
    """

    def check_dot_oam_9dots_helper(E):
        flagPlus, flagMinus = True, True
        minIndex = np.argmin(E)
        for i in range(minIndex - len(E), minIndex - 1, 1):
            if E[i] >= E[i + 1]:
                flagMinus = False
                break
        maxIndex = np.argmax(E)
        for i in range(maxIndex - len(E), maxIndex - 1, 1):
            if E[i] <= E[i + 1]:
                flagPlus = False
                break
        if flagPlus:
            # print(np.arg() + np.arg() - np.arg() - np.arg())
            return True, +1
        elif flagMinus:
            return True, -1
        return False, 0

    shape = np.shape(E)
    ans = np.zeros(shape)
    for i in range(1, shape[0] - 1, 1):
        for j in range(1, shape[1] - 1, 1):
            Echeck = np.array([E[i - 1, j - 1], E[i - 1, j], E[i - 1, j + 1],
                               E[i, j + 1], E[i + 1, j + 1], E[i + 1, j],
                               E[i + 1, j - 1], E[i, j - 1]])
            oamFlag, oamValue = check_dot_oam_9dots_helper(Echeck)
            if oamFlag:
                ######
                ans[i - circle:i + 1 + circle, j - circle:j + 1 + circle] = nonValue
                #####
                ans[i, j] = oamValue * value
                if bigSingularity:
                    ans[i - 1:i + 2, j - 1:j + 2] = oamValue * value
            else:
                ans[i, j] = nonValue
    return ans


def plane_singularities_finder_4dots(E, circle, value, nonValue, bigSingularity):
    """
    cut_non_oam helper. see cut_non_oam for more details
    """

    def check_dot_oam_4dots_helper(E):
        def arg(x):
            return np.angle(np.exp(1j * x))

        sum = arg(E[1] - E[0]) + arg(E[2] - E[3]) - arg(E[2] - E[1]) - arg(E[1] - E[0])
        if sum > 3:
            return True, +1
        if sum < -3:
            return True, -1
        return False, 0

    shape = np.shape(E)
    ans = np.zeros(shape)
    for i in range(1, shape[0] - 1, 1):
        for j in range(1, shape[1] - 1, 1):
            Echeck = np.array([E[i, j], E[i, j + 1], E[i + 1, j + 1], E[i + 1, j]])
            oamFlag, oamValue = check_dot_oam_4dots_helper(Echeck)
            if oamFlag:
                ######
                ans[i - circle:i + 1 + circle, j - circle:j + 1 + circle] = nonValue
                #####
                ans[i, j] = oamValue * value
                if bigSingularity:
                    ans[i - 1:i + 2, j - 1:j + 2] = oamValue * value
            else:
                ans[i, j] = nonValue
    return ans


def fill_dict_as_matrix_helper(E, dots=None, nonValue=0, check=False):
    """
    cut_non_oam helper. see cut_non_oam for more details
    """
    if dots is None:
        dots = {}
    shape = np.shape(E)
    if len(shape) == 3:
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    if E[i, j, k] != nonValue:
                        if check:
                            if dots.get((i, j, k)) is None:
                                dots[(i, j, k)] = E[i, j, k]
                        else:
                            dots[(i, j, k)] = E[i, j, k]
    else:
        for i in range(shape[0]):
            for j in range(shape[1]):
                if E[i, j] != nonValue:
                    if check:
                        if dots.get((i, j, 0)) is None:
                            dots[(i, j, 0)] = E[i, j]
                    else:
                        dots[(i, j, 0)] = E[i, j]
    return dots


def cut_non_oam(E, value=1, nonValue=0, bigSingularity=False, axesAll=False, circle=1,
                singularities_finder=plane_singularities_finder_9dots):
    """
    this function finds singularities
    returns [3D Array, dots only]
    :param E: complex field
    :param value: all singularities will have these values +- values (depend on the sing sign)
    :param nonValue: all non-singularities have this value
    :param bigSingularity: singularities and all dots around it has "value"
    :param axesAll: singularities are searched not only in Oxy, but in Oxz and Oyz
    :param circle: if the singularity is found, the circle around it is automatically equaled to nonValue
    :param singularities_finder: plane_singularities_finder_9dots or _4dots. 2nd one is faster
    :return: [the whole 3d massive with values and nonValues, dict[x, y, z]=value]
    """
    shape = np.shape(E)
    if len(shape) == 2:
        ans = singularities_finder(E, circle, value, nonValue, bigSingularity)
        ans[:1, :] = nonValue
        ans[-1:, :] = nonValue
        ans[:, :1] = nonValue
        ans[:, -1:] = nonValue
        dots = fill_dict_as_matrix_helper(ans)
    else:
        ans = np.copy(E)
        for i in range(shape[2]):
            ans[:, :, i] = cut_non_oam(ans[:, :, i], value=value, nonValue=nonValue,
                                       bigSingularity=bigSingularity)[0]
        dots = fill_dict_as_matrix_helper(ans)

        if axesAll:
            for i in range(shape[1]):
                ans[:, i, :] += cut_non_oam(E[:, i, :], value=value, nonValue=nonValue,
                                            bigSingularity=bigSingularity)[0]
            dots = fill_dict_as_matrix_helper(ans, dots, check=True)
            for i in range(shape[0]):
                ans[i, :, :] += cut_non_oam(E[i, :, :], value=value, nonValue=nonValue,
                                            bigSingularity=bigSingularity)[0]
            dots = fill_dict_as_matrix_helper(ans, dots, check=True)
    # print(ans)
    return ans, dots


def get_singularities(E, value=1, nonValue=0, bigSingularity=False, axesAll=False, circle=1,
                      singularities_finder=plane_singularities_finder_4dots):
    """
    cut_non_oam simplifier. Just return the array of singularities
    """
    if isinstance(E, dict):
        dotsOnly = E
    else:
        dotsFull, dotsOnly = cut_non_oam(E, value, nonValue, bigSingularity, axesAll, circle,
                                         singularities_finder)
    dots = np.array([list(dots) for (dots, OAM) in dotsOnly.items()])
    return dots


def W_energy(EArray, xArray=None, yArray=None):
    """
    total power in Oxy plane
    :param EArray:
    :param xArray:
    :param yArray:
    :return:
    """
    if xArray is None or yArray is None:
        shape = np.shape(EArray)
        xArray = np.arange(shape[0])
        yArray = np.arange(shape[1])
    dx = xArray[1] - xArray[0]
    dy = yArray[1] - yArray[0]
    W = np.real(np.sum(np.conj(EArray) * EArray) * dx * dy)
    return W


def Jz_calc_no_conj(EArray, xArray=None, yArray=None):
    EArray = np.array(EArray)
    Er, Ei = np.real(EArray), np.imag(EArray)
    if xArray is None or yArray is None:
        shape = np.shape(EArray)
        xArray = np.arange(shape[0])
        yArray = np.arange(shape[1])
    x0 = (xArray[-1] + xArray[0]) / 2
    y0 = (yArray[-1] + yArray[0]) / 2
    x = np.array(xArray) - x0
    y = np.array(yArray) - y0
    dx = xArray[1] - xArray[0]
    dy = yArray[1] - yArray[0]
    sumJz = 0
    for i in range(1, len(xArray) - 1, 1):
        for j in range(1, len(yArray) - 1, 1):
            dErx = (Er[i + 1, j] - Er[i - 1, j]) / (2 * dx)
            dEry = (Er[i, j + 1] - Er[i, j - 1]) / (2 * dy)
            dEix = (Ei[i + 1, j] - Ei[i - 1, j]) / (2 * dx)
            dEiy = (Ei[i, j + 1] - Ei[i, j - 1]) / (2 * dy)
            # dErx = (Er[i + 1, j] - Er[i, j]) / (dx)
            # dEry = (Er[i, j + 1] - Er[i, j]) / (dy)
            # dEix = (Ei[i + 1, j] - Ei[i, j]) / (dx * 2)
            # dEiy = (Ei[i, j + 1] - Ei[i, j]) / (dy)
            # print(x[i] * Er[i, j] * dEiy, - y[j] * Er[i, j] * dEix, -
            #           x[i] * Ei[i, j] * dEry, + y[j] * Ei[i, j] * dErx)
            sumJz += (x[i] * Er[i, j] * dEiy - y[j] * Er[i, j] * dEix -
                      x[i] * Ei[i, j] * dEry + y[j] * Ei[i, j] * dErx)
    # Total moment
    Jz = (sumJz * dx * dy)
    W = W_energy(EArray)
    print(f'Total OAM charge = {Jz / W}\tW={W}')
    return Jz


def integral_number2_OAMcoefficients_FengLiPaper(fieldFunction, r, l):
    """
    Implementation of the Integral (2) from the FengLi paper for calculating the weight of OAM in r
    :param fieldFunction:
    :param r: radius where you want to know OAM
    :param l: exp(1j * j * phi)
    """

    # function helps to get y value from x and r. Sign is used in 2 different parts of the CS.
    # helper => it is used only in other functions, you don't use it yourself
    def y_helper(x, sign, r):
        return sign * np.sqrt(r ** 2 - x ** 2)

    def f1(x):  # function f in the upper half - plane
        Y = y_helper(x, +1, r)
        return fieldFunction(x, Y) * (-1) / Y * np.exp(-1j * l * fg.phi(x, Y)) / np.sqrt(2 * np.pi)

    def f2(x):  # function f in the lower half - plane
        Y = y_helper(x, -1, r)
        return fieldFunction(x, Y) * (-1) / Y * np.exp(-1j * l * fg.phi(x, Y)) / np.sqrt(2 * np.pi)

    i1 = fg.integral_of_function_1D(f1, r, -r)  # upper half - plane integration
    i2 = fg.integral_of_function_1D(f2, -r, r)  # lower half - plane integration
    answer = i1[0] + i2[0]  # [0] - is the integral value, [1:] - errors value and other stuff, we don't need it
    return answer


def integral_number3_OAMpower_FengLiPaper(fieldFunction, rMin, rMax, rResolution, l):
    """
    Implementation of the Integral (3) from the FengLi paper
    Calculating total power in the OAM with charge l
    :param fieldFunction: function of the field
    :param rMin, rMax: boundaries of the integral
    :param l: OAM
    :return: Pl
    """
    rArray = np.linspace(rMin, rMax, rResolution)
    aRArray = np.zeros(rResolution, dtype=complex)
    for i in range(rResolution):
        aRArray[i] = integral_number2_OAMcoefficients_FengLiPaper(fieldFunction, rArray[i], l)
    pL = integrate.simps(np.abs(aRArray) ** 2 * rArray, rArray)  # using interpolation
    return pL


class Singularities3D:
    """
    Work with singularities of any 3D complex field
    """

    def __init__(self, field3D=None):
        """
        :param field3D: any 3D complex field
        """
        self.field3D = field3D
        self.dotsDict = None  # self.dotsXY or self.dotsAll (can switch with self.swap()
        self.dotsXY = None  # singularities from XY planes
        self.dotsAll = None  # singularities from XY+XZ+YZ planes
        self.dotsList = None  # np.array [[x,y,z], [x,y,z], ...] random order
        self.mesh = None  # np.meshgrid from LG_combination
        self.coefficients = None  # [Cl1p1, Cl2p2...] from LG_combination
        self.modes = None  # [(l1,p1), (l2,p2) ...] from LG_combination
        # self.fill_dotsDict_from_field3D(_dotsXY=True)

    def field_LG_combination(self, mesh, coefficients, modes, **kwargs):
        """
        creating the field of any combination of LG beams
        Sum(Cl1p1 * LG_simple(*mesh, l=l1, p=p1, **kwargs))
        :param mesh: np.meshgrid
        :param coefficients: [Cl1p1, Cl2p2...] ...
        :param modes: [(l1,p1), (l2,p2) ...]
        """
        field = 0
        self.mesh = mesh
        self.coefficients = coefficients
        self.modes = modes
        for num, coefficient in enumerate(coefficients):
            field += coefficient * bp.LG_simple(*mesh, l=modes[num][0], p=modes[num][1], **kwargs)
        self.field3D = field

    def plot_plane_2D(self, zPlane, show=True, **kwargs):
        """
        Plot any z plane, both abs and angle
        :param zPlane: number of the plane to plot (0<=z<=shape[2])
        :return: None
        """
        pl.plot_2D(np.abs(self.field3D[:, :, zPlane]), **kwargs)
        pl.plot_2D(np.angle(self.field3D[:, :, zPlane]), map='hsv', **kwargs)
        if show:
            plt.show()

    def plot_center_2D(self, **kwargs):
        """
        Plot the center plane (z=0 if from z is from -l to l)
        :return: None
        """
        shape = np.shape(self.field3D)
        self.plot_plane_2D(shape[2] // 2, **kwargs)

    def plot_dots(self, show=True, **kwargs):
        """
        Plot self.dots (scatter) using fOAM.plot_knot_dots()
        if self.dots is not initialized, initialization with self.fill_dotsDict_from_field3D()
        :param kwargs: Everything for fOAM.plot_knot_dots()
         also for self.fill_dotsDict_from_field3D()
        :return: None
        """
        if self.dotsDict is None:
            self.fill_dotsDict_from_field3D(**kwargs)
        ax = plot_knot_dots(self.dotsDict, **kwargs)
        if show:
            plt.show()
        return ax
        # fg.distance_between_points()

    def plot_density(self, **kwargs):
        """
        Plot density on the browser
        :kwargs: Everything for fOAM.plot_knot_dots()
        :return: None
        """
        pl.plot_3D_density(np.angle(self.field3D), **kwargs)

    def fill_dotsDict_from_field3D(self, _dotsXY=True, **kwargs):
        """
        Filing self.dots with self.dotsXY. for self.dotsALL use parameter _dotsXY
        :param kwargs: Everything for fg.cut_non_oam()
        :param _dotsXY: if True, we are filling with self.dotsXY, otherwise with self.dotsALL
        :return: number of dots in self.dots
        """
        if _dotsXY:
            if self.dotsXY is None:
                self.fill_dotsXY(**kwargs)
            self.dotsDict = self.dotsXY
        else:
            if self.dotsAll is None:
                self.fill_dotsAll(**kwargs)
            self.dotsDict = self.dotsAll
        return len(self.dotsDict)

    def fill_dotsList(self):
        self.dotsList = np.array([[x, y, z] for (x, y, z) in self.dotsDict])

    def fill_dotsXY(self, **kwargs):
        """
        fill in self.dotsXY with using only XY cross-sections for singularities
        :param kwargs: fg.cut_non_oam besides axesAll
        :return:
        """
        garbage, self.dotsXY = cut_non_oam(np.angle(self.field3D), axesAll=False, **kwargs)
        self.dotsDict = self.dotsXY

    def fill_dotsAll(self, **kwargs):
        """
        fill in self.dotsALL with using ALL 3 cross-sections for singularities
        :param kwargs: fg.cut_non_oam besides axesAll
        :return:
        """
        garbage, self.dotsAll = cut_non_oam(np.angle(self.field3D), axesAll=True, **kwargs)
        self.dotsDict = self.dotsAll

    def dots_swap(self, **kwargs):
        """
        change self.dots between self.dotsXY and self.dotsAll
        if self.dots is not either of those -> self.fill_dotsDict_from_field3D()
        print the new self.dots
        :return: None
        """
        if self.dotsDict is self.dotsXY:
            if self.dotsAll is None:
                self.fill_dotsAll()
            self.dotsDict = self.dotsAll
            print(f'Dots are now in all 3 planes')
        elif self.dotsDict is self.dotsAll:
            if self.dotsXY is None:
                self.fill_dotsXY()
            self.dotsDict = self.dotsXY
            print(f'Dots are now in the XY-plane')
        else:
            self.fill_dotsDict_from_field3D(*kwargs)
            print(f'Dots were not dotsXY or dotsAll. Now dots are in the XY-plane')


class Knot(Singularities3D):
    """
    Knot field (unknots also are knots in that case)
    """

    def __init__(self, field3D=None):
        """
        :param field3D: any 3D complex field
        """
        Singularities3D.__init__(self, field3D)
        self.dotsKnotList = None  # the actual knot (ordered line)
        self.knotSP = None

    def build_knot_pyknotid(self, **kwargs):
        """
        function build normilized pyknotid knot
        :return:
        """
        if self.dotsKnotList is None:
            self.fill_dotsKnotList()
        zMid = (max(z for x, y, z in self.dotsKnotList) + min(z for x, y, z in self.dotsKnotList)) / 2
        xMid = (max(x for x, y, z in self.dotsKnotList) + min(x for x, y, z in self.dotsKnotList)) / 2
        yMid = (max(y for x, y, z in self.dotsKnotList) + min(y for x, y, z in self.dotsKnotList)) / 2
        self.knotSP = sp.Knot(np.array(self.dotsKnotList) - [xMid, yMid, zMid], add_closure=False, **kwargs)

    def plot_knot(self, **kwargs):
        """
        plot the knot
        """
        if self.dotsKnotList is None:
            self.fill_dotsKnotList()
        if self.knotSP is None:
            self.build_knot_pyknotid(**kwargs)
        plt.plot([1], [1])
        self.knotSP.plot()
        plt.show()

    def fill_dotsKnotList(self):
        if self.dotsList is None:
            self.fill_dotsList()
        distance_matrix = euclidean_distance_matrix(self.dotsList, self.dotsList)
        permutation, distance = solve_tsp_local_search(distance_matrix)
        # print(dots[permutation])
        # print(permutation)
        self.dotsKnotList = self.dotsList[permutation]

    def fill_dotsKnotList_mine(self):
        """
        fill in self.dotsList by removing charge sign and placing everything into the list [[x, y, z], [x, y, z]...]
        :return: None
        """

        def min_dist(dot, dots):
            elements = [(fg.distance_between_points(dot, d), i) for i, d in enumerate(dots)]
            minEl = min(elements, key=lambda i: i[0])
            return minEl

        self.dotsKnotList = []
        dotsDict = {}
        for [x, y, z] in self.dotsDict:
            if not (z in dotsDict):
                dotsDict[z] = []
            dotsDict[z].append([x, y])
        indZ = next(iter(dotsDict))  # z coordinate
        indInZ = 0  # dot number in XY plane at z
        indexes = np.array([-1, 0, 1])  # which layers we are looking at
        currentDot = dotsDict[indZ].pop(indInZ)
        # distCheck = 20
        while dotsDict:
            # print(indZ, currentDot, dotsDict)
            minList = []  # [min, layer, position in Layer] for all indexes + indZ layers
            for i in indexes + indZ:  # searching the closest element among indexes + indZ
                if not (i in dotsDict):
                    continue
                minVal, min1Ind = min_dist(currentDot, dotsDict[i])
                # if minVal <= distCheck:
                minList.append([minVal, i, min1Ind])
            if not minList:
                newPlane = 2
                while not minList:
                    for i in [-newPlane, newPlane] + indZ:  # searching the closest element among indexes + indZ
                        if not (i in dotsDict):
                            continue
                        minVal, min1Ind = min_dist(currentDot, dotsDict[i])
                        # if minVal <= distCheck:
                        minList.append([minVal, i, min1Ind])
                    newPlane += 1
                if newPlane > 3:
                    print(f'we have some dots left. Stopped')
                    print(indZ, currentDot, dotsDict)
                    break
                print(f'dots are still there, the knot builred cannot use them all\nnew plane: {newPlane}')
            minFin = min(minList, key=lambda i: i[0])
            # if minFin[1] != indZ:
            self.dotsKnotList.append([*dotsDict[minFin[1]].pop(minFin[2]), minFin[1]])
            currentDot = self.dotsKnotList[-1][:-1]  # changing the current dot to a new one
            indZ = minFin[1]
            # else:
            #     dotsDict[minFin[1]].pop(minFin[2])
            #     currentDot = self.dotsList[-1][:-1]  # changing the current dot to a new one
            #     indZ = minFin[1]
            # currentDot = self.dotsList[-1][:-1][:]  # changing the current dot to a new one
            # indZ = minFin[1]
            # dotsDict[minFin[1]].pop(minFin[2])
            if not dotsDict[indZ]:  # removing the empty plane (0 dots left)
                del dotsDict[indZ]

    def check_knot_alex(self) -> bool:
        checkVal = None
        if self.knotSP is None:
            self.build_knot_pyknotid()
        t = sympy.symbols("t")
        self.alexPol = self.knotSP.alexander_polynomial(variable=t)
        if self.__class__.__name__ == 'Trefoil':
            checkVal = -t ** 2 + t - 1
        if checkVal is None:
            print(f'There is no check value for this type of knots')
            return False
        if self.alexPol == checkVal:
            return True
        return False


if __name__ == '__main__':
    import timeit
    import sympy
    from python_tsp.distances import euclidean_distance_matrix
    from python_tsp.heuristics import solve_tsp_local_search


    def func_time_main():
        xMinMax = 3
        yMinMax = 3
        zMinMax = 0.8
        zRes = 70
        xRes = yRes = 80
        xyzMesh = fg.create_mesh_XYZ(xMinMax, yMinMax, zMinMax, xRes, yRes, zRes, zMin=None)
        coeff = [1.715, -5.662, 6.381, -2.305, -4.356]
        phase = [0, 0, 0, 0, 0]
        coeff = [a * np.exp(1j * p) for a, p in zip(coeff, phase)]
        beam = bp.LG_combination(*xyzMesh, coefficients=coeff, modes=((0, 0), (0, 1), (0, 2), (0, 3), (3, 0)))
        # dots = get_singularities(np.angle(beam), bigSingularity=False, axesAll=False)
        # pl.plot_2D(np.abs(beam[:,:, zRes//2]))
        dotsExp = np.load('C:\\Users\\Dima\\Box\\Knots Exp\\'
                          'Experimental Data\\dots\\trefoil\\Field SR = 0.95\\3foil_turb_25.npy',
                          allow_pickle=True).item()
        dots = get_singularities(dotsExp)
        # pl.plot_scatter_3D(dots[:, 0], dots[:, 1], dots[:, 2])
        distance_matrix = euclidean_distance_matrix(dots, dots)
        minSum = 0
        checkValue = 3
        checkNumber = 3
        indStay = []
        for i, line in enumerate(distance_matrix):
            value = line[line > 0].min()
            minSum += value
            lineSorted = np.sort(line)
            if lineSorted[checkNumber] < checkValue:
                indStay.append(i)
            # print(i, value)
        # print(indStay, print(len(indStay)))
        # print(minSum / np.shape(distance_matrix)[0])
        newDots = dots[indStay]
        pl.plot_scatter_3D(newDots[:, 0], newDots[:, 1], newDots[:, 2])
        for i, line in enumerate(distance_matrix):
            value = line[line > 0].min()
            minSum += value
            if value < checkValue:
                indStay.append(i)
            print(i, value)
        exit()
        permutation, distance = solve_tsp_local_search(distance_matrix)
        print(permutation[:10])
        print(distance[:10])
        # print(dots[permutation])
        # print(permutation)
        dotsKnotList = dots[permutation]
        exit()
        # plot_knot_dots(beam, show=True)


    timeit.timeit(func_time_main, number=1)
