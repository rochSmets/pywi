
import sys
import os

sys.path.append(os.path.join(os.environ['PYWI_DIR'], 'runs'))

import numpy as np
import heckle


def saddles(mat : np.ndarray) -> list:

    """
    returns the list of all saddle points of the input matrix

    """


    (N, M) = mat.shape

    jMax = np.argmax(mat, axis = 1) # index of col for max in each row
    iMin = np.argmin(mat, axis = 0) # index of row for min in each col

    IJMax = [(i,jMax[i]) for i in range(N)] # list of indexes of max of each row
    IJMin = [(iMin[j],j) for j in range(M)] # list of indexes of min of each col

    maxRowMinCol = list(set(IJMax) & set(IJMin)) # max of row, min of col


    iMax = np.argmax(mat, axis = 0) # index of row for max in each col
    jMin = np.argmin(mat, axis = 1) # index of col for min in each row

    IJMax = [(iMax[j],j) for j in range(M)] # list of indexes of max of each col
    IJMin = [(i,jMin[i]) for i in range(N)] # list of indexes of min of each row

    minRowMaxCol = list(set(IJMax) & set(IJMin)) # min of row, max of col


    return maxRowMinCol + minRowMaxCol





def reconnectedFlux(path: str, time : list, xPointBox: list, axis : int) -> list:

    """

    """

    run = heckle.Heckle('/run/media/smets/croq0Disk/blAckDog/02b/', 'Nb')

    # index of the box where to find X point
    ijBox = np.array(np.array(xPointBox)/run.dl, dtype=int)

    # index of guessed X point (from xPointBox)
    ijG = np.mean(ijBox, axis=1, dtype=int)

    times, groups = run.getTimeGroups(time)
    times.sort()

    norm = lambda x, y : (y[0]-x[0])^2+(y[1]-x[1])^2

    timeOut = []
    ijXOut = []
    fluxOut = []
    EzOut = []

    for index, time in enumerate(times) :
        print('time : {0:6.2f} / {1:6.2f}'.format(time, times[-1]))

        Az = run.fourierFlux(time)
        az = Az[ijBox[0][0]:ijBox[0][1],ijBox[1][0]:ijBox[1][1]]

        saddlePoints = saddles(az)

        minDist = np.inf

        for i in range(len(saddlePoints)):
            # print("X point # {0:1} : [{1:4} - {2:4}]".format(i, saddlePoints[i][0]+ijBox[0][0], saddlePoints[i][1]+ijBox[1][0]))
            ijS = [saddlePoints[i][0]+ijBox[0][0], saddlePoints[i][1]+ijBox[1][0]]
            ijS = ijG ##########################################################
            dist = norm(ijG, ijS)

            # keep the best X point, ie closest to the guessed one
            if dist < minDist:
                ijX = ijS
                minDist = dist

        AzLim = 0.5*(Az[ijX[0], 0]+Az[ijX[0], -1])
        flux = np.fabs(Az[ijX[0], ijX[1]]-AzLim)
        Ez = run.GetE(time)[ijX[0], ijX[1] ,2]

        # print(time, ijX, flux, Ez)

        timeOut.append(time)
        ijXOut.append(ijX)
        fluxOut.append(flux)
        EzOut.append(Ez)

    return timeOut, ijXOut, fluxOut, EzOut

