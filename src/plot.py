# plot.py
from __future__ import division

import travelTime as tt
import numpy as np
import matplotlib.pyplot as plt


def cleanDisplay():
    plt.close()
    plt.figure(figsize=(10, 7))


def showDisplay():
    plt.title("")
    plt.xlabel("Time [ns]", fontsize=16, labelpad=16)
    plt.ylabel("Reflection coefficient [-]", fontsize=16, labelpad=16)
    plt.tick_params(axis='both', which='major', labelsize=12, pad=8)
    plt.tick_params(axis='both', which='minor', labelsize=12, pad=8)
    plt.ylim(bottom=-1, top=1)
    plt.show()


def drawWaveForm(showDerivatives):
    lastIndex = len(tt.reflecCoeff) - 2
    t = np.zeros(lastIndex, float)
    for i in range(lastIndex):
        t[i] = tt.timeVector[i] * 1E09
    y = tt.reflecCoeff[0:lastIndex]
    dy = tt.dy[0:lastIndex]
    dy2 = tt.dy2[0:lastIndex]
    plt.plot(t, y, 'k.')
    if showDerivatives:
        plt.plot(t, dy, 'k--')
        plt.plot(t, dy2, 'r--')


def drawRegressionLines():
    nrPoints = len(tt.timeVector)
    step = int(16. * (nrPoints / 256.0))

    t = np.zeros(nrPoints, float)
    line1 = np.zeros(nrPoints, float)
    line2 = np.zeros(nrPoints, float)
    flatLine = np.zeros(nrPoints, float)
    for i in range(nrPoints):
        t[i] = tt.timeVector[i] * 1E09
        line1[i] = tt.line1.a * tt.timeVector[i] + tt.line1.b
        line2[i] = tt.line2.a * tt.timeVector[i] + tt.line2.b
        flatLine[i] = tt.lastFlatLine.b

    index = int(tt.p2.x / tt.deltaTime)
    first = max(0, index - step)
    last = min(nrPoints, index + step)
    plt.plot(t[first:last], line1[first:last], 'k')
    plt.plot(t[first:last], line2[first:last], 'k')

    first = nrPoints - step * 4
    last = nrPoints
    plt.plot(t[first:last], flatLine[first:last], 'k')

    plt.plot(tt.p0.x * 1E09, tt.p0.y, 'rs')
    plt.plot(tt.p1.x * 1E09, tt.p1.y, 'rs')
    plt.plot(tt.p2.x * 1E09, tt.p2.y, 'rs')
    plt.plot(tt.p3.x * 1E09, tt.p3.y, 'rs')
