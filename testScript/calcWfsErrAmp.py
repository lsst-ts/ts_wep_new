import os

from lsst.ts.wep.WEPController import WEPController, plotDonutImg

from lsst.ts.wep.SciWFDataCollector import SciWFDataCollector
from lsst.ts.wep.SciIsrWrapper import SciIsrWrapper
from lsst.ts.wep.SourceSelector import SourceSelector
from lsst.ts.wep.SourceProcessor import SourceProcessor, abbrevDectectorName
from lsst.ts.wep.WFEstimator import WFEstimator
from lsst.ts.wep.Utility import getModulePath

from cwfs.Tool import plotImage

if __name__ == "__main__":

    # This script is to calculate the wavefront error from the amplifier images.
    # DM team only supports the LSST FAM mode and SE only supports the ComCam and 
    # WFS at this moment.
    # Therefore, in this script, the WFS data is LSST central raft and the parameters 
    # of calculating the wavefront error is in the ComCam 1.5 mm condition.
    
    # Instintiate the components
    sourSelc = SourceSelector()
    dataCollector = SciWFDataCollector()
    isrWrapper = SciIsrWrapper()
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
    homeDir = os.path.expanduser("~")
    pathOfRawData = os.path.join(homeDir, "Document", "phosimObsData", "raw")
    destinationPath = os.path.join(homeDir, "Document", "phosimObsData", "input")
    dataCollector.config(pathOfRawData=pathOfRawData, destinationPath=destinationPath)

    # Configurate the ISR wrapper
    postIsrImgDir = os.path.join(homeDir, "Document", "phosimObsData", "output")
    isrWrapper.configWrapper(inputs=destinationPath, outputs=postIsrImgDir)
    isrWrapper.configBulter(postIsrImgDir)

    # Configurate the source processor
    focalPlaneFolder = os.path.join(modulePath, "test")
    sourProc.config(donutRadiusInPixel=starRadiusInPixel, folderPath2FocalPlane=focalPlaneFolder, 
                    pixel2Arcsec=0.2)

    # Configurate the wavefront estimator
    defocalDisInMm = 1.5
    
    # Size of donut in pixel if defocal distance = 1.5 mm
    sizeInPix = 160
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
    skyInfoFilePath = os.path.join(homeDir, "Document", "phosimObsData", "skyInfo", "skyLsstFamInfo.txt")

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
    # The data butler can only import the unique data id for one time
    for obsId in obsIdList:
        fitsFileArg = "lsst_*%d*.fits.gz" % obsId
        try:
            wepCntlr.ingestSimImages(fitsFileArg=fitsFileArg)
        except RuntimeError:
            print("Data is registered in the data butler already.")

    # Do the ISR
    for obsId in obsIdList:
        for sensorName in sensorNameList:
            wepCntlr.doISR(obsId, sensorName)

    # Get the wfs images
    # It is noted that the image of assembled CCD is in DM coordinate.
    # Therefor, set "expInDmCoor=True" to rotate to camera team coordinate.
    wfsImgMap = wepCntlr.getPostISRDefocalImgMap(sensorNameList, obsIdList=obsIdList, 
                                                 expInDmCoor=True)

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