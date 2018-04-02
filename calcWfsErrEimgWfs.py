import os

from wep.WEPController import WEPController, plotDonutImg
from wep.WFDataCollector import WFDataCollector
from wep.SourceSelector import SourceSelector
from wep.SourceProcessor import SourceProcessor, abbrevDectectorName
from wep.WFEstimator import WFEstimator

from cwfs.Tool import plotImage

if __name__ == "__main__":
    
    # This script is to calculate the wavefront error from the eimages of corner WFS.
    # This helps the evaluation of wavefront error without the silicon model.
    
    # Instintiate the components
    sourSelc = SourceSelector()
    dataCollector = WFDataCollector()
    sourProc = SourceProcessor()

    instruFolderPath = os.path.join(".", "instruData")
    algoFolderPath = os.path.join(".", "algo")
    wfsEsti = WFEstimator(instruFolderPath, algoFolderPath)

    # Configurate the source selector
    cameraType = "lsst"
    dbType = "LocalDb"
    aFilter = "g"
    cameraMJD = 59580.0

    sourSelc.configSelector(cameraType=cameraType, dbType=dbType, aFilter=aFilter, 
                            cameraMJD=cameraMJD)

    # Set the criteria of neighboring stars
    starRadiusInPixel = 63
    spacingCoefficient = 2.5
    sourSelc.configNbrCriteria(starRadiusInPixel, spacingCoefficient)

    # Configurate the WFS data collector
    # Data butler does not support the corner WFS at this moment.
    pathOfRawData = os.path.join(".", "test", "phosimOutput")
    destinationPath = butlerInputs = butlerOutputs = os.path.join(".", "test")
    dataCollector.config(pathOfRawData=pathOfRawData, destinationPath=destinationPath)

    # Configurate the source processor
    focalPlaneFolder = os.path.join(".", "test")
    sourProc.config(donutRadiusInPixel=starRadiusInPixel, folderPath2FocalPlane=focalPlaneFolder, 
                    pixel2Arcsec=0.2)

    # Configurate the wavefront estimator
    defocalDisInMm = None
    
    # Size of donut in pixel for corner WFS
    sizeInPix = 120
    wfsEsti.config(solver="exp", instName=cameraType, opticalModel="offAxis", 
                    defocalDisInMm=defocalDisInMm, sizeInPix=sizeInPix)

    # Initiate the WEP Controller
    wepCntlr = WEPController()
    wepCntlr.config(sourSelc=sourSelc, dataCollector=dataCollector, sourProc=sourProc, 
                    wfsEsti=wfsEsti)

    # Set the database address
    dbAdress = os.path.join(".", "test", "bsc.db3")

    # Do the query
    pointing = (0,0)
    cameraRotation = 0.0
    skyInfoFilePath = os.path.join(".", "test", "phosimOutput", "realWfs", "output", 
                                    "skyWfsInfo.txt")

    camOrientation = "corner"
    neighborStarMap, starMap, wavefrontSensors = wepCntlr.getTargetStarByFile(dbAdress, 
                                            skyInfoFilePath, pointing, cameraRotation, 
                                            orientation=camOrientation, tableName="TempTable")

    # Get the available sensor name list
    sensorNameList = list(starMap.keys())

    # Import the PhoSim simulated image
    # Because the data butler does not support the corner WFS at this moment,
    # the steps related to it such as doISR() have been ignored.
    wfsDir = os.path.join("realWfs", "output")
    cornerWfsImgMap = wepCntlr.getPostISRDefocalImgMap(sensorNameList, wfsDir=wfsDir)
    wfsImgMap = wepCntlr.getPostISRDefocalImgMap(sensorNameList, wfsDir=wfsDir)

    # Get the donut images
    donutMap = wepCntlr.getDonutMap(neighborStarMap, wfsImgMap, aFilter, doDeblending=False, 
                                    sglDonutOnly=True)

    # Plot the donut images
    saveToDir = os.path.join(".", "test", "donutImg")
    plotDonutImg(donutMap, saveToDir=saveToDir, dpi=None)

    # Calculate the wavefront error for the individual donut
    donutMap = wepCntlr.calcWfErr(donutMap)
    for sensorName, donutImgList in donutMap.items():
        for donutImg in donutImgList:
            if (donutImg.zer4UpNm is not None):
                print(sensorName)
                print(donutImg.starId)
                print(donutImg.zer4UpNm)

    # Calculate the average wavefront error
    for sensorName, donutImgList in donutMap.items():
        if (not sensorName.endswith("B")):
            avgErr = wepCntlr.calcSglAvgWfErr(donutImgList)

            print(sensorName)
            print(avgErr)