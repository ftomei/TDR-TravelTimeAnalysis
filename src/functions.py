# functions.py

from __future__ import division
from math import sqrt, exp

NODATA = -9999
SPEED_OF_LIGHT = 299792458  # [m s^-1]
airPermittivity = 1.00058986  # [-]
waterDensity = 997  # [kg m^-3]


class WaveForm:
    probeLenght = NODATA
    travelTime = NODATA
    vp = NODATA
    v1 = NODATA
    vf = NODATA
    vr = NODATA
    wcTopp = NODATA
    wcMalicky = NODATA
    wcMixModel = NODATA
    wcCurioni = NODATA


def getLiquidPermittivity(temperature):
    deltaT = temperature - 25.
    return 78.54 * (1 - 4.579E-03 * deltaT)


def getBulkPermittivity(waveForm):
    return ((SPEED_OF_LIGHT * waveForm.vp * waveForm.travelTime) / (2. * waveForm.probeLenght)) ** 2


def getWaterContentTopp(bulkPermittivity):
    return -5.3E-02 + 2.92E-02 * bulkPermittivity - 5.5E-04 * bulkPermittivity ** 2 + 4.3E-06 * bulkPermittivity ** 3


def getWaterContentMalicki(bulkPermittivity, bulkDensity):
    bulkDensity /= 1000.
    return (sqrt(bulkPermittivity) - 0.819 - 0.168 * bulkDensity - 0.159 * bulkDensity ** 2) / (
                7.17 + 1.18 * bulkDensity)


def getWaterContentMixModel(bulkPermittivity, bulkDensity,
                            solidPermittivity, liquidPermittivity, alpha):
    porosity = 1. - bulkDensity / 2650.
    numerator = bulkPermittivity ** alpha - ((1. - porosity) * solidPermittivity ** alpha
                                             + porosity * airPermittivity ** alpha)
    denominator = liquidPermittivity ** alpha - airPermittivity ** alpha
    return numerator / denominator


def getBulkDensityJung(bulkPermittivity, v1, vf, c1, d1, f1):
    numerator = v1 / vf
    denominator = c1 + d1 * (bulkPermittivity - 1) - c1 * exp(-f1 * (bulkPermittivity - 1))
    return waterDensity * (numerator / denominator)


def getBulkDensityCurioni(waveForm, a, b, c):
    bulkPermittivity = getBulkPermittivity(waveForm)
    numerator = waterDensity * waveForm.vr
    denominator = a + b * (waveForm.v1 * sqrt(bulkPermittivity)) ** c
    return numerator / denominator


def getWaterContentCurioni(bulkPermittivity, bulkDensity, a, b):
    return (1 / b) * (sqrt(bulkPermittivity) * (waterDensity / bulkDensity) - a) / 100.
