import os

from lsst.ts.wep.WEPController import WEPController, plotDonutImg

from lsst.ts.wep.WFDataCollector import WFDataCollector
from lsst.ts.wep.EimgIsrWrapper import EimgIsrWrapper
from lsst.ts.wep.SourceSelector import SourceSelector
from lsst.ts.wep.SourceProcessor import SourceProcessor, abbrevDectectorName
from lsst.ts.wep.WFEstimator import WFEstimator
from lsst.ts.wep.Utility import getModulePath

from cwfs.Tool import plotImage

if __name__ == "__main__":
    
    # This script is to calculate the wavefront error from the eimages of ComCam.
    # This helps the evaluation of wavefront error without the silicon model.
    
    # Instintiate the components
    sourSelc = SourceSelector()
    dataCollector = WFDataCollector()
    isrWrapper = EimgIsrWrapper()
    sourProc = SourceProcessor()

    modulePath = getModulePath()
    instruFolderPath = os.path.join(modulePath, "algoData", "cwfs", "instruData")
    algoFolderPath = os.path.join(modulePath, "algoData", "cwfs", "algo")
    wfsEsti = WFEstimator(instruFolderPath, algoFolderPath)

    # Configurate the source selector
    cameraType = "comcam"
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
    pathOfRawData = os.path.join(modulePath, "test", "phosimOutput")
    destinationPath = butlerInputs = butlerOutputs = os.path.join(modulePath, "test")
    regisAdress = os.path.join(modulePath, "test", "registry.sqlite3")
    dataCollector.config(pathOfRawData=pathOfRawData, destinationPath=destinationPath, 
                dbAdress=regisAdress, butlerInputs=butlerInputs, butlerOutputs=butlerOutputs)

    # Configurate the ISR wrapper
    isrWrapper.configWrapper(inputs=butlerInputs, outputs=butlerOutputs)
    isrWrapper.configBulter(butlerInputs, outputs=butlerOutputs)

    # Configurate the source processor
    focalPlaneFolder = os.path.join(modulePath, "test")
    sourProc.config(donutRadiusInPixel=starRadiusInPixel, folderPath2FocalPlane=focalPlaneFolder, 
                    pixel2Arcsec=0.2)

    # Configurate the wavefront estimator
    defocalDisInMm = 1
    
    # Size of donut in pixel if defocal distance = 1 mm
    sizeInPix = 120
    wfsEsti.config(solver="exp", instName=cameraType, opticalModel="offAxis", 
                    defocalDisInMm=defocalDisInMm, sizeInPix=sizeInPix)

    # Initiate the WEP Controller
    wepCntlr = WEPController()
    wepCntlr.config(sourSelc=sourSelc, dataCollector=dataCollector, isrWrapper=isrWrapper, 
                    sourProc=sourProc, wfsEsti=wfsEsti)

    # Set the database address
    dbAdress = os.path.join(modulePath, "test", "bsc.db3")

    # Do the query
    pointing = (0,0)
    cameraRotation = 0.0
    skyInfoFilePath = os.path.join(modulePath, "test", "phosimOutput", "realComCam2", "output", 
                                   "skyComCamInfo.txt")

    camOrientation = "all"
    neighborStarMap, starMap, wavefrontSensors = wepCntlr.getTargetStarByFile(dbAdress, 
                                            skyInfoFilePath, pointing, cameraRotation, 
                                            orientation=camOrientation, tableName="TempTable")

    # Get the available sensor name list
    sensorNameList = list(starMap.keys())

    # Observation IDs
    extraObsId = 9005000
    intraObsId = 9005001
    obsIdList = [intraObsId, extraObsId]

    # Import the PhoSim simulated image
    intraImgDir = os.path.join("realComCam2", "output", "Intra")
    extraImgDir = os.path.join("realComCam2", "output", "Extra")
    dataDirList = [intraImgDir, extraImgDir]
    for dataDir in dataDirList:
        wepCntlr.ingestSimImages(dataDir=dataDir, atype="raw", overwrite=False)

    # Do the ISR
    for obsId in obsIdList:
        for sensorName in sensorNameList:
            wepCntlr.doISR(obsId, sensorName)

    # Get the wfs images
    # It is noted that the eimage is in camera team coordinate
    wfsImgMap = wepCntlr.getPostISRDefocalImgMap(sensorNameList, obsIdList=obsIdList, 
                                                 expInDmCoor=False)

    # Get the donut images
    donutMap = wepCntlr.getDonutMap(neighborStarMap, wfsImgMap, aFilter, doDeblending=False, 
                                    sglDonutOnly=True)

    # Plot the donut images
    saveToDir = os.path.join(modulePath, "test", "donutImg")
    plotDonutImg(donutMap, saveToDir=saveToDir, dpi=None)

    # Calculate the wavefront error for the individual donut
    donutMap = wepCntlr.calcWfErr(donutMap)
    for sensorName, donutImgList in donutMap.items():
        for donutImg in donutImgList:
            if (donutImg.zer4UpNm is not None):
                print(sensorName)
                print(donutImg.starId)
                print(donutImg.zer4UpNm)

    # Generate the master donut images
    masterDonutMap = wepCntlr.generateMasterImg(donutMap)

    # Plot the master donut image
    for sensorName, img in masterDonutMap.items():
        fileName = abbrevDectectorName(sensorName)

        intraFilePath = os.path.join(saveToDir, fileName + "_intra.png")
        plotImage(img[0].intraImg, title=fileName+"_intra", show=False, saveFilePath=intraFilePath)
        
        extraFilePath = os.path.join(saveToDir, fileName + "_extra.png")
        plotImage(img[0].extraImg, title=fileName+"_extra", show=False, saveFilePath=extraFilePath)

    # Calculate the wavefront error for the master donut
    masterDonutMap = wepCntlr.calcWfErr(masterDonutMap)
    for sensorName, masterDonutImg in masterDonutMap.items():
        print(sensorName)
        print(masterDonutImg[0].zer4UpNm)

    # Calculate the average wavefront error
    for sensorName, donutImgList in donutMap.items():
        if (not sensorName.endswith("B")):
            avgErr = wepCntlr.calcSglAvgWfErr(donutImgList)

            print(sensorName)
            print(avgErr)