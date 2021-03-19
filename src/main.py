# ------------------------------------------------------------------
# TDR TRAVEL TIME
# ------------------------------------------------------------------
# 2020.10 
# Authors:
# Fausto Tomei      ftomei@arpae.it
# Marco Bittelli    marco.bittelli@unibo.it
# ------------------------------------------------------------------
from __future__ import division

from readDataFile import readDataFile, readGenericDataFile
from functions import *
from Marquardt import Marquardt
from plot import *

import os
import sys

if sys.version_info >= (3, 0):
    from tkinter import *  # 3.x
    from tkinter.filedialog import askopenfilename, askdirectory
    from tkinter.messagebox import showerror
else:
    from Tkinter import *  # 2.7
    from tkFileDialog import askopenfilename, askdirectory
    from tkMessageBox import showerror

isDataLoaded = False
waterTemperature = 20
liquidPermittivity = getLiquidPermittivity(waterTemperature)

mainWindow = Tk()
mainWindow.title("TDRPy")
mainWindow.geometry("%dx%d" % (500, 600))

headerNrStr = StringVar()
vpStr = StringVar()
probeLenghtStr = StringVar()
windowBeginStr = StringVar()
windowWidthStr = StringVar()
probeHandleStr = StringVar()
handlePermittivityStr = StringVar()
point0XStr = StringVar()
point1XStr = StringVar()
point2XStr = StringVar()
v1Str = StringVar()
vfStr = StringVar()
ratioStr = StringVar()
estBulkDensityStr = StringVar()
ttStr = StringVar()
bulkPermittivityStr = StringVar()
waterTemperatureStr = StringVar()
bulkDensityStr = StringVar()
aStr = StringVar()
bStr = StringVar()
cStr = StringVar()
solidPermittivityStr = StringVar()
liquidPermittivityStr = StringVar()
geometricParStr = StringVar()
wcToppStr = StringVar()
wcMalickiStr = StringVar()
wcMixModelStr = StringVar()
wcCurioniStr = StringVar()

if sys.version_info >= (3, 0):
    waterTempLabel = Label(mainWindow, text="Water temp. [\u00B0C]")
else:
    waterTempLabel = Label(mainWindow, text="Water temp. [C]")
waterTempLabel.place(x=50, y=295)
waterTempWidget = Entry(mainWindow, width=6, textvariable=waterTemperatureStr)
waterTempWidget.place(x=190, y=295)


def getEpsilonLabel():
    if sys.version_info >= (3, 0):
        return Label(mainWindow, font="helvetica 12", text="\u03F5")
    else:
        return Label(mainWindow, font="helvetica 10", text="e")


def checkTWater(event):
    global waterTemperature, liquidPermittivity
    if event.widget == waterTempWidget:
        waterTemperature = float(waterTemperatureStr.get())
        liquidPermittivity = getLiquidPermittivity(waterTemperature)
        liquidPermittivityStr.set(format(liquidPermittivity, ".3f"))


def importData():
    global isDataLoaded
    nrHeaderValues = int(headerNrStr.get())
    fileName = askopenfilename()
    if fileName != "":
        x, isFileOk = readDataFile(fileName, nrHeaderValues, '\t', False)
        if not isFileOk:
            showerror("Wrong file", "Wrong file!\nError reading row nr." + str(x))
            return False
        if len(x) == 1:
            data = x[0, :]
        else:
            data = x[:, 0]
        tt.reflecCoeff = tt.normalizeVector(data)
        print("number of values:", len(data))
        isDataLoaded = True
        ComputeTT(True)


def ComputeTT(isShow=True):
    if not isDataLoaded:
        showerror("Warning", "Data are not loaded.")
        return False

    # read parameters
    vp = float(vpStr.get())
    # fraction of speed of light [-]
    probeLenght = float(probeLenghtStr.get())
    windowBegin = float(windowBeginStr.get())
    windowWidth = float(windowWidthStr.get())
    probeHandle = float(probeHandleStr.get())
    handlePermittivity = float(handlePermittivityStr.get())
    bulkDensity = float(bulkDensityStr.get())
    solidPermittivity = float(solidPermittivityStr.get())
    geomParameter = float(geometricParStr.get())
    a = float(aStr.get())
    b = float(bStr.get())
    c = float(cStr.get())

    # compute
    nrPoints = len(tt.reflecCoeff)
    tt.WF_parameters(vp, probeHandle, windowBegin, windowWidth, nrPoints)

    if not tt.computeTravelTime(probeHandle, handlePermittivity, vp):
        showerror("Warning", "Wrong data, header or parameter")
        return

    w = WaveForm()
    w.travelTime = tt.p2.x - tt.p1.x
    w.probeLenght = probeLenght
    w.vp = vp
    bulkPermittivity = getBulkPermittivity(w)
    w.wcTopp = getWaterContentTopp(bulkPermittivity)
    w.wcMalicki = getWaterContentMalicki(bulkPermittivity, bulkDensity)
    w.wcMixModel = getWaterContentMixModel(bulkPermittivity, bulkDensity, solidPermittivity,
                                           liquidPermittivity, geomParameter)

    w.v1 = abs(tt.p2.y - tt.p1.y)
    w.vf = abs(tt.p3.y + 1)
    w.vr = w.v1 / w.vf
    w.bdCurioni = getBulkDensityCurioni(w, a, b, c)
    w.wcCurioni = getWaterContentCurioni(bulkPermittivity, w.bdCurioni, a, b)

    # print results
    x0 = tt.p0.x * 1E09
    point0XStr.set(format(x0, '.3f'))
    x1 = tt.p1.x * 1E09
    point1XStr.set(format(x1, '.3f'))
    x2 = tt.p2.x * 1E09
    point2XStr.set(format(x2, '.3f'))

    ttStr.set(format(w.travelTime * 1E09, '.3f'))
    bulkPermittivityStr.set(format(bulkPermittivity, '.2f'))

    v1Str.set(format(w.v1, '.3f'))
    vfStr.set(format(w.vf, '.3f'))
    ratioStr.set(format(w.vr, '.3f'))
    estBulkDensityStr.set(format(w.bdCurioni, ".1f"))

    wcToppStr.set(format(w.wcTopp, ".3f"))
    wcMalickiStr.set(format(w.wcMalicki, ".3f"))
    wcMixModelStr.set(format(w.wcMixModel, ".3f"))
    wcCurioniStr.set(format(w.wcCurioni, ".3f"))

    if isShow:
        cleanDisplay()
        drawWaveForm(False)
        drawRegressionLines()
        showDisplay()

    return w


def computeAllWaveform():
    global isDataLoaded

    myPath = askdirectory() + "/"
    nrHeaderValues = int(headerNrStr.get())

    # output file
    outFile = open(myPath + "output.csv", "w")
    outFile.write("soil code, bulk permittivity, v1, vr, bulk density, "
                  "wc Topp, wc Malicki, wc MixModel, wc Curioni\n")

    # cycle on files
    files = os.listdir(myPath)
    for fileName in files:
        suffix = fileName[len(fileName) - 3:]
        if suffix == "dat":
            soilCode = fileName[:len(fileName) - 4]
            fileName = myPath + fileName
            x, isFileOk = readDataFile(fileName, nrHeaderValues, '\t', False)
            if isFileOk:
                if len(x) == 1:
                    data = x[0, :]
                else:
                    data = x[:, 0]
                tt.reflecCoeff = tt.normalizeVector(data)
                isDataLoaded = True
                w = ComputeTT(False)

                # compute water content with Curioni bulk density
                # CHECK
                bulkPermittivity = getBulkPermittivity(w)
                w.wcTopp = getWaterContentTopp(bulkPermittivity)
                w.wcMalicki = getWaterContentMalicki(bulkPermittivity, w.bdCurioni)

                solidPermittivity = float(solidPermittivityStr.get())
                geomParameter = float(geometricParStr.get())
                w.wcMixModel = getWaterContentMixModel(bulkPermittivity, w.bdCurioni, solidPermittivity,
                                                       liquidPermittivity, geomParameter)
                # output
                outFile.write("%s,%.2f,%.3f,%.3f,%.1f,%.3f,%.3f,%.3f,%.3f\n" % (
                    soilCode, bulkPermittivity, w.v1, w.vr, w.bdCurioni, w.wcTopp,
                    w.wcMalicki, w.wcMixModel, w.wcCurioni))

    outFile.close()


def calibrateCurioni():
    global isDataLoaded
    myPath = askdirectory() + "/"
    obsFileName = myPath + "obs_density.csv"
    if not os.path.exists(obsFileName):
        showerror("Missing file", "It requires observed density (obs_density.csv)")
        return False

    # Read obs data
    x, isFileOk = readGenericDataFile(obsFileName, 1, ',', False)
    fileList = []
    obsDensityAll = np.zeros(len(x))
    soilListAll = []
    for i in range(len(x)):
        soilListAll.append(x[i][0])
        fileList.append(soilListAll[i] + ".dat")
        obsDensityAll[i] = x[i][1]

    # read waveform
    waveForm = []
    soilList = []
    obsDensity = []
    nrHeaderValues = int(headerNrStr.get())
    files = os.listdir(myPath)
    for fileName in files:
        if fileName in fileList:
            i = fileList.index(fileName)
            fileName = myPath + fileName
            x, isFileOk = readDataFile(fileName, nrHeaderValues, '\t', False)
            if isFileOk:
                if len(x) == 1:
                    data = x[0, :]
                else:
                    data = x[:, 0]
                tt.reflecCoeff = tt.normalizeVector(data)
                isDataLoaded = True

                w = ComputeTT(False)
                waveForm.append(w)
                soilList.append(soilListAll[i])
                obsDensity.append(obsDensityAll[i])

    # Marquardt
    a = float(aStr.get())
    b = float(bStr.get())
    c = float(cStr.get())
    v0 = np.array([a, b, c], float)
    vmin = np.array([0, 0, 0], float)
    vmax = np.array([0.1, 1.0, 10.0], float)
    v, estDensity = Marquardt(v0, vmin, vmax, waveForm, obsDensity)
    aStr.set("{0:.4f}".format(v[0]))
    bStr.set("{0:.4f}".format(v[1]))
    cStr.set("{0:.4f}".format(v[2]))

    # output
    outFile = open(myPath + "fitting.csv", "w")
    outFile.write("%.4f,%.4f,%.4f\n" % (v[0], v[1], v[2]))
    outFile.write("soil code, bulk permittivity, v1, vr, obs density, est density\n")
    for i in range(len(obsDensity)):
        k = getBulkPermittivity(waveForm[i])
        outFile.write("%s,%.2f,%.3f,%.3f,%.1f,%.1f\n" % (
            soilList[i], k, waveForm[i].v1, waveForm[i].vr, obsDensity[i], estDensity[i]))
    outFile.close()


def main():
    vpStr.set(0.99)
    probeLenghtStr.set(0.1)
    windowBeginStr.set(0.)
    windowWidthStr.set(3.)
    probeHandleStr.set(0.120)
    handlePermittivityStr.set(1.7)
    point0XStr.set(0)
    point1XStr.set(0)
    point2XStr.set(0)

    v1Str.set(0)
    vfStr.set(0)
    ratioStr.set(0)
    estBulkDensityStr.set(0)
    ttStr.set(0)
    bulkPermittivityStr.set(0)
    bulkDensityStr.set(1350)
    waterTemperatureStr.set(waterTemperature)
    solidPermittivityStr.set(4.0)
    liquidPermittivityStr.set(format(liquidPermittivity, ".3f"))
    geometricParStr.set(0.5)
    aStr.set(0.0)
    bStr.set(0.14)
    cStr.set(1.2)
    wcToppStr.set(0)
    wcMalickiStr.set(0)
    wcMixModelStr.set(0)
    wcCurioniStr.set(0)

    posX_col1 = 50
    posY_col1 = 80
    posX_col2 = posX_col1 + 240
    posY_col2 = posY_col1

    buttonImport = Button(mainWindow, text="Import data", command=importData)
    buttonImport.place(x=posX_col1, y=15)

    headerLabel = Label(mainWindow, text="Header values nr.")
    headerLabel.place(x=posX_col1 + 80, y=20)
    headerWidget = Entry(mainWindow, width=3, textvariable=headerNrStr)
    headerWidget.place(x=posX_col1 + 160, y=20)
    headerWidget.insert(0, "8")

    computeTTButton = Button(mainWindow, text="  Compute  ", command=ComputeTT)
    computeTTButton.place(x=posX_col1 + 220, y=15)

    vpLabel = Label(mainWindow, text="Vp [-]")
    vpLabel.place(x=posX_col1, y=posY_col1)
    vpWidget = Entry(mainWindow, width=6, textvariable=vpStr)
    vpWidget.place(x=posX_col1 + 140, y=posY_col1)
    posY_col1 += 25

    probeLenghtLabel = Label(mainWindow, text="Probe length [m]")
    probeLenghtLabel.place(x=posX_col1, y=posY_col1)
    probeLenghtWidget = Entry(mainWindow, width=6, textvariable=probeLenghtStr)
    probeLenghtWidget.place(x=posX_col1 + 140, y=posY_col1)
    posY_col1 += 25

    winBeginLabel = Label(mainWindow, text="Window begin [m]")
    winBeginLabel.place(x=posX_col1, y=posY_col1)
    winBeginWidget = Entry(mainWindow, width=6, textvariable=windowBeginStr)
    winBeginWidget.place(x=posX_col1 + 140, y=posY_col1)
    posY_col1 += 25

    winwidthLabel = Label(mainWindow, text="Window width [m]")
    winwidthLabel.place(x=posX_col1, y=posY_col1)
    winwidthWidget = Entry(mainWindow, width=6, textvariable=windowWidthStr)
    winwidthWidget.place(x=posX_col1 + 140, y=posY_col1)
    posY_col1 += 25

    probeHandleLabel = Label(mainWindow, text="Probe handle [m]")
    probeHandleLabel.place(x=posX_col1, y=posY_col1)
    probeHandleWidget = Entry(mainWindow, width=6, textvariable=probeHandleStr)
    probeHandleWidget.place(x=posX_col1 + 140, y=posY_col1)
    posY_col1 += 25

    epsilonLabel = getEpsilonLabel()
    epsilonLabel.place(x=posX_col1, y=posY_col1-6)
    handleLabel = Label(mainWindow, text="handle")
    handleLabel.place(x=posX_col1+10, y=posY_col1)
    permittivityWidget = Entry(mainWindow, width=6, textvariable=handlePermittivityStr)
    permittivityWidget.place(x=posX_col1 + 140, y=posY_col1)
    posY_col1 += 40

    # Soil parameters
    SoilParameterLabel = Label(mainWindow, text="Soil parameters", font="helvetica 10 bold", fg="red")
    SoilParameterLabel.place(x=posX_col1 - 20, y=posY_col1)
    posY_col1 += 25

    bulkDensityLabel = Label(mainWindow, text="Bulk density [kg m-3]")
    bulkDensityLabel.place(x=posX_col1, y=posY_col1)
    bulkDensityWidget = Entry(mainWindow, width=6, textvariable=bulkDensityStr)
    bulkDensityWidget.place(x=posX_col1 + 140, y=posY_col1)
    posY_col1 += 50

    epsilon2Label = getEpsilonLabel()
    epsilon2Label.place(x=posX_col1, y=posY_col1-5)
    epsilon2Label = Label(mainWindow, text="liquid")
    epsilon2Label.place(x=posX_col1 + 10, y=posY_col1)
    liquidPermittivityLabel = Label(mainWindow, textvariable=liquidPermittivityStr)
    liquidPermittivityLabel.place(x=posX_col1 + 140, y=posY_col1)
    posY_col1 += 25

    epsilon3Label = getEpsilonLabel()
    epsilon3Label.place(x=posX_col1, y=posY_col1-5)
    epsilon3Label = Label(mainWindow, text="solid")
    epsilon3Label.place(x=posX_col1 + 10, y=posY_col1)
    solidPermittivityWidget = Entry(mainWindow, width=6, textvariable=solidPermittivityStr)
    solidPermittivityWidget.place(x=posX_col1 + 140, y=posY_col1)
    posY_col1 += 25

    geometricParLabel = Label(mainWindow, text="alpha (geom. param.)")
    geometricParLabel.place(x=posX_col1, y=posY_col1)
    geometricParWidget = Entry(mainWindow, width=6, textvariable=geometricParStr)
    geometricParWidget.place(x=posX_col1 + 140, y=posY_col1)
    posY_col1 += 40

    # travel time output 
    ttResultsLabel = Label(mainWindow, text="Travel Time results", font="helvetica 10 bold", fg="blue")
    ttResultsLabel.place(x=posX_col1 - 20, y=posY_col1)
    posY_col1 += 25

    point0Label = Label(mainWindow, text="point 0 x [ns]")
    point0Label.place(x=posX_col1, y=posY_col1)
    point0Widget = Label(mainWindow, width=6, textvariable=point0XStr)
    point0Widget.place(x=posX_col1 + 140, y=posY_col1)
    posY_col1 += 25

    point1Label = Label(mainWindow, text="point 1 x [ns]")
    point1Label.place(x=posX_col1, y=posY_col1)
    point1Widget = Label(mainWindow, width=6, textvariable=point1XStr)
    point1Widget.place(x=posX_col1 + 140, y=posY_col1)
    posY_col1 += 25

    point2Label = Label(mainWindow, text="point 2 x [ns]")
    point2Label.place(x=posX_col1, y=posY_col1)
    point2Widget = Label(mainWindow, width=6, textvariable=point2XStr)
    point2Widget.place(x=posX_col1 + 140, y=posY_col1)
    posY_col1 += 25

    ttLabel = Label(mainWindow, text="Travel Time [ns]")
    ttLabel.place(x=posX_col1, y=posY_col1)
    ttWidget = Label(mainWindow, width=6, textvariable=ttStr)
    ttWidget.place(x=posX_col1 + 140, y=posY_col1)
    posY_col1 += 25

    bulkPermittivityLabel = Label(mainWindow, text="Bulk permittivity")
    bulkPermittivityLabel.place(x=posX_col1, y=posY_col1)
    bulkPermittivityWidget = Label(mainWindow, width=6, textvariable=bulkPermittivityStr)
    bulkPermittivityWidget.place(x=posX_col1 + 140, y=posY_col1)

    # Column 2: Curioni parameter
    bdParamLabel = Label(mainWindow, text="Curioni parameters", font="helvetica 10 bold", fg="red")
    bdParamLabel.place(x=posX_col2 - 20, y=posY_col2)
    buttonCurioni = Button(mainWindow, text="Fitting", command=calibrateCurioni)
    buttonCurioni.place(x=posX_col2 + 120, y=posY_col2 - 2)
    posY_col2 += 35

    aLabel = Label(mainWindow, text="a [-]")
    aLabel.place(x=posX_col2, y=posY_col2)
    aWidget = Entry(mainWindow, width=6, textvariable=aStr)
    aWidget.place(x=posX_col2 + 120, y=posY_col2)
    posY_col2 += 25

    bLabel = Label(mainWindow, text="b [-]")
    bLabel.place(x=posX_col2, y=posY_col2)
    bWidget = Entry(mainWindow, width=6, textvariable=bStr)
    bWidget.place(x=posX_col2 + 120, y=posY_col2)
    posY_col2 += 25

    cLabel = Label(mainWindow, text="c [-]")
    cLabel.place(x=posX_col2, y=posY_col2)
    cWidget = Entry(mainWindow, width=6, textvariable=cStr)
    cWidget.place(x=posX_col2 + 120, y=posY_col2)
    posY_col2 += 30

    computeAll = Button(mainWindow, text="Compute all waveform", command=computeAllWaveform)
    computeAll.place(x=posX_col2, y=posY_col2)
    posY_col2 += 50

    # Bulk density
    densityLabel = Label(mainWindow, text="Bulk density", font="helvetica 10 bold", fg="blue")
    densityLabel.place(x=posX_col2 - 20, y=posY_col2)
    posY_col2 += 25

    v1Label = Label(mainWindow, text="V1 [-]")
    v1Label.place(x=posX_col2, y=posY_col2)
    v1Widget = Label(mainWindow, width=6, textvariable=v1Str)
    v1Widget.place(x=posX_col2 + 120, y=posY_col2)
    posY_col2 += 25

    vfLabel = Label(mainWindow, text="Vf [-]")
    vfLabel.place(x=posX_col2, y=posY_col2)
    vfWidget = Label(mainWindow, width=6, textvariable=vfStr)
    vfWidget.place(x=posX_col2 + 120, y=posY_col2)
    posY_col2 += 25

    ratioLabel = Label(mainWindow, text="V1 / Vf [-]")
    ratioLabel.place(x=posX_col2, y=posY_col2)
    ratioWidget = Label(mainWindow, width=6, textvariable=ratioStr)
    ratioWidget.place(x=posX_col2 + 120, y=posY_col2)
    posY_col2 += 25

    estBulkDensityLabel = Label(mainWindow, text="Bulk density [kg m3]")
    estBulkDensityLabel.place(x=posX_col2, y=posY_col2)
    estBulkDensityWidget = Label(mainWindow, width=6, textvariable=estBulkDensityStr)
    estBulkDensityWidget.place(x=posX_col2 + 120, y=posY_col2)
    posY_col2 += 25
    posY_col2 += 40

    # Water content
    wcLabel = Label(mainWindow, text="Water Content", font="helvetica 10 bold", fg="blue")
    wcLabel.place(x=posX_col2 - 20, y=posY_col2)
    posY_col2 += 25

    ToppLabel = Label(mainWindow, text="Topp")
    ToppLabel.place(x=posX_col2, y=posY_col2)
    ToppWidget = Label(mainWindow, width=6, textvariable=wcToppStr)
    ToppWidget.place(x=posX_col2 + 120, y=posY_col2)
    posY_col2 += 25

    MalickiLabel = Label(mainWindow, text="Malicki")
    MalickiLabel.place(x=posX_col2, y=posY_col2)
    MalickiWidget = Label(mainWindow, width=6, textvariable=wcMalickiStr)
    MalickiWidget.place(x=posX_col2 + 120, y=posY_col2)
    posY_col2 += 25

    dielecMixModelLabel = Label(mainWindow, text="Diel. mix model")
    dielecMixModelLabel.place(x=posX_col2, y=posY_col2)
    dielecMixModelWidget = Label(mainWindow, width=6, textvariable=wcMixModelStr)
    dielecMixModelWidget.place(x=posX_col2 + 120, y=posY_col2)
    posY_col2 += 25

    CurioniLabel = Label(mainWindow, text="Curioni")
    CurioniLabel.place(x=posX_col2, y=posY_col2)
    CurioniWcWidget = Label(mainWindow, width=6, textvariable=wcCurioniStr)
    CurioniWcWidget.place(x=posX_col2 + 120, y=posY_col2)

    mainWindow.bind("<Leave>", checkTWater)
    mainWindow.protocol("WM_DELETE_WINDOW", mainWindow.destroy)

    mainWindow.mainloop()


main()
