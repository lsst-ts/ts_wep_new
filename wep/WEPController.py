import os, re
import numpy as np

import matplotlib
# Must be before importing matplotlib.pyplot or pylab!
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from astropy.io import fits
from scipy.ndimage.measurements import center_of_mass

from wep.WFDataCollector import WFDataCollector
from wep.IsrWrapper import getImageData
from wep.EimgIsrWrapper import EimgIsrWrapper
from wep.SourceSelector import SourceSelector
from wep.SourceProcessor import SourceProcessor, abbrevDectectorName
from wep.WFEstimator import WFEstimator

from wep.LocalDatabaseDecorator import LocalDatabaseDecorator
from wep.DefocalImage import DefocalImage, DonutImage

class WEPController(object):

    def __init__(self):
        
        self.sourSelc = None
        self.dataCollector = None
        self.isrWrapper = None
        self.sourProc = None
        self.wfsEsti = None
        self.middleWare = None

    def config(self, sourProc=None, dataCollector=None, isrWrapper=None, sourSelc=None, 
                wfsEsti=None, middleWare=None):

        self.__setVar(sourSelc, "sourSelc")        
        self.__setVar(dataCollector, "dataCollector")
        self.__setVar(isrWrapper, "isrWrapper")
        self.__setVar(sourProc, "sourProc")
        self.__setVar(wfsEsti, "wfsEsti")
        self.__setVar(middleWare, "middleWare")

    def getWfsList(self):
        """
        
        Get the corner wavefront sensor (WFS) list in the canonical form.
        
        Returns:
            [list] -- WFS name list.
        """

        wfsList = ["R:0,0 S:2,2,A", "R:0,0 S:2,2,B", "R:0,4 S:2,0,A", "R:0,4 S:2,0,B", 
                   "R:4,0 S:0,2,A", "R:4,0 S:0,2,B", "R:4,4 S:0,0,A", "R:4,4 S:0,0,B"]

        return wfsList

    def __setVar(self, value, attrName):
        """
        
        Set the value of attribute.
        
        Arguments:
            value {[obj]} -- New value.
            attrName {[str]} -- Attribute name to set the value.
        """

        if (value is not None):
            setattr(self, attrName, value)

    def getTargetStarByFile(self, dbAdress, skyInfoFilePath, pointing, cameraRotation, 
                            orientation=None, tableName="TempTable"):
        """
        
        Get the target stars by querying the file.
        
        Arguments:
            dbAdress {[str]} -- Local database address.
            skyInfoFilePath {[str]} -- File path of sky information.
            pointing {[tuple]} -- Camera boresight (RA, Decl) in degree.
            cameraRotation {[float]} -- Camera rotation angle in degree.
        
        Keyword Arguments:
            orientation {[str]} -- Orientation of wavefront sensor(s) on camera. (default: {None})
            tableName {[str]} -- Table name. (default: {None})
        
        Returns:
            {[dict]} -- Information of neighboring stars and candidate stars with the name of 
                        sensor as a dictionary.
            {[dict]} -- Information of stars with the name of sensor as a dictionary.
            {[dict]} -- Corners of sensor with the name of sensor as a dictionary.
        """

        # Check the database name is local database
        if (self.sourSelc.name != self.sourSelc.LocalDb):
            raise TypeError("The database type is not LocalDatabaseDecorator.")

        # Get the filter type
        aFilter = self.sourSelc.getFilter()

        # Connect the database
        self.sourSelc.connect(dbAdress)

        # Create the table
        self.sourSelc.db.createTable(aFilter, tableName)

        # Insert the sky data
        self.sourSelc.db.insertDataByFile(aFilter, tableName, skyInfoFilePath, skiprows=1)
        
        # Do the query and analysis
        neighborStarMap, starMap, wavefrontSensors = self.sourSelc.getTargetStar(pointing, cameraRotation, 
                                                                orientation=orientation, tableName=tableName)

        neighborStarMap, starMap, wavefrontSensors = self.__analyzeStarMap(neighborStarMap, starMap, 
                                                                                    wavefrontSensors)

        # Delete the table
        self.sourSelc.db.deleteTable(tableName)
        
        # Disconnect the database
        self.sourSelc.disconnect()

        return neighborStarMap, starMap, wavefrontSensors

    def __analyzeStarMap(self, neighborStarMap, starMap, wavefrontSensors):
        """
        
        Analyze the star map and remove the sensor without bright stars.
        
        Arguments:
            neighborStarMap {[dict]} -- Information of neighboring stars and candidate stars with 
                                        the name of sensor as a dictionary.
            starMap {[dict]} -- Information of stars with the name of sensor as a dictionary.
            wavefrontSensors {[dict]} -- Corners of sensor with the name of sensor as a dictionary.
        
        Returns:
            {[dict]} -- Information of neighboring stars and candidate stars with the name 
                        of sensor as a dictionary.
            {[dict]} -- Information of stars with the name of sensor as a dictionary.
            {[dict]} -- Corners of sensor with the name of sensor as a dictionary.
        """

        # Collect the sensor list without the bright star
        noStarSensorList = []
        for aKey, aItem in starMap.items():
            if len(aItem.RA) == 0:
                noStarSensorList.append(aKey)

        # Remove the keys in map
        for aKey in noStarSensorList:
            neighborStarMap.pop(aKey)
            starMap.pop(aKey)
            wavefrontSensors.pop(aKey)
        
        return neighborStarMap, starMap, wavefrontSensors

    def importPhoSimDataToButler(self, dataDir, atype="raw", overwrite=False):
        """
        
        Import the PhoSim simulated data to match with the data butler to use. This means the 
        registry.sqlite3 repo will be inserted with the meta data if necessary.
        
        Arguments:
            dataDir {[str]} -- PhoSim FITS data directory.
        
        Keyword Arguments:
            atype {[str]} -- Dataset type. (default: {"raw"})
            overwrite {[boolean]} -- Overwrite the existed files or not. (default: {False})
        
        Raises:
            ValueError -- Not allowed type ("raw", "bias", "dark", "flat").
        """

        # Get all files in the directory
        fileList = self.__getRawFileList(dataDir)

        # Find the obsId and aFilter
        obsIdList = []
        for fileName in fileList:
            m = re.match(r"\S*_(\d*)_f(\d)_\S*", fileName)
            if (m is not None):
                obsIdList.append(m.groups())

        # Get the unique list
        obsIdList = list(set(obsIdList))
        if (len(obsIdList) > 1):
            raise RuntimeError("There are more than one unique ObdId and filter in directory.")
        data = obsIdList[0]

        # Import to butler
        phosimFilterID = {"0": "u", "1": "g", "2": "r", "3": "i", "4": "z", "5": "y"}
        obsId = int(data[0])
        aFilter = phosimFilterID[data[1]]
        self.dataCollector.importPhoSimDataToButler(dataDir, obsId=obsId, aFilter=aFilter, atype=atype, 
                                                    overwrite=overwrite)

    def __getRawFileList(self, dataDir):
        """
        
        Get the raw file list in the directory.
        
        Arguments:
            dataDir {[str]} -- Data directory.
        
        Returns:
            [list] -- File list.
        """

        # Get all files in the directory
        fullDataDir = os.path.join(self.dataCollector.pathOfRawData, dataDir)
        fileList = [f for f in os.listdir(fullDataDir) if os.path.isfile(os.path.join(fullDataDir, f))]

        return fileList

    def doISR(self, visit, sensorName, snap=0, fakeDatasetType="eimage", 
                outputDatasetType="postISRCCD"):
        """
        
        Do the instrument signature removal (ISR).
        
        Arguments:
            visit {[int]} -- Visit time.
            sensorName {[str]} -- Sensor name. (e.g. "R:2,2 S:1,1")
        
        Keyword Arguments:
            snap {int} -- Snap time (0 or 1) means first/ second exposure. (default: {0})
            fakeDatasetType {[str]} -- Use this type of image supported by lsst camera mapper 
                                        to simulate the post-ISR image. (default: {"eimage"})
            outputDatasetType {[str]} -- Output data type supported by lsst camera mapper. 
                                        (default: {"postISRCCD"})

        Returns:
            [ExposureU] -- Exposure image after ISR.
        """

        # Use the regular expression to analyze the input name
        raft, sensor, channel = self.__getSensorInfo(sensorName)
        if (raft is not None):

            # Do the ISR
            self.isrWrapper.doISR(visit, snap, raft, sensor, channel=None, 
                        fakeDatasetType=fakeDatasetType, outputDatasetType=outputDatasetType)
        else:
            raise RuntimeError("Sensor name: '%s' is not allowed." % sensorName)


    def __getSensorInfo(self, sensorName):
        """
        
        Get the sensor information.
        
        Arguments:
            sensorName {[str]} -- Sensor name (e.g. "R:2,2 S:1,1" or "R:0,0 S:2,2,A")
        
        Returns:
            [str] -- Raft.
            [str] -- Sensor.
            [str] -- Channel.
        """

        raft = sensor = channel = None
        
        # Use the regular expression to analyze the input name
        m = re.match(r"R:(\d,\d) S:(\d,\d)(?:,([A,B]))?$", sensorName)
        if (m is not None):
            raft, sensor, channel = m.groups()[0:3]

        return raft, sensor, channel

    def __searchFileName(self, fileList, matchName, snap=0):
        """
        
        Search the file name in list.
        
        Arguments:
            fileList {[list]} -- File name list.
            matchName {[str]} -- Match name.

        Keyword Arguments:
            snap {int} -- Snap number (default: {0})
        
        Returns:
            [str] -- Matched file name.
        """

        matchFileName = None
        for fileName in fileList:
            m = re.match(r"\S*%s_E00%d\S*" % (matchName, snap), fileName)

            if (m is not None):
                matchFileName = m.group()
                break

        return matchFileName

    def getPostISRDefocalImgMap(self, sensorNameList, obsIdList=None, wfsDir=None, snap=0):
        """
        
        Get the post-ISR defocal image map.
        
        Arguments:
            sensorNameList {[list]} -- List of sensor name which is in the canonical form.
        
        Keyword Arguments:
            obsIdList {[list]} -- Observation Id list in [intraObsId, extraObsId]. (default: {None})
            wfsDir {[str]} -- Directory to wavefront sensor image data. (default: {None})
            snap {int} -- Snap number (default: {0})
        
        Returns:
            [dict] -- Post-ISR image map.
        """

        # Construct the dictionary 
        wfsImgMap = {}

        # Get the file list
        if (wfsDir is not None):
            # Get the file list
            wfsFileList = self.__getRawFileList(wfsDir)

        # Get the waveront image map
        for sensorName in sensorNameList:

            # Get the sensor name information
            raft, sensor, channel = self.__getSensorInfo(sensorName)

            # The intra/ extra defocal images are decided by obsId
            if (obsIdList is not None):

                imgList = []
                for ii in range(2):
                    dataId = dict(visit=obsIdList[ii], snap=snap, raft=raft, sensor=sensor)
                    img = self.isrWrapper.butler.get(datasetType="postISRCCD", dataId=dataId, 
                                                     immediate=True)
                    imgList.append(getImageData(img))

                wfsImgMap[sensorName] = DefocalImage(intraImg=imgList[0], extraImg=imgList[1])

            # The intra/ extra defocal images are decided by physical configuration
            # C0: intra, C1: extra
            if (wfsDir is not None):
                
                # Get the abbreviated name
                abbrevName = abbrevDectectorName(sensorName)

                # Search for the file name
                matchFileName = self.__searchFileName(wfsFileList, abbrevName, snap=snap)
                
                if (matchFileName is not None):
                    
                    # Get the file name
                    fitsFilsPath = os.path.join(self.dataCollector.pathOfRawData, wfsDir, 
                                                matchFileName)
                    wfsImg = fits.getdata(fitsFilsPath)

                    # Add image to map
                    wfsImgMap[sensorName] = DefocalImage()

                    # "C0" = "A" = "Intra-focal image"
                    if (channel=="A"):
                        wfsImgMap[sensorName].setImg(intraImg=wfsImg)
                    # "C1" = "B" = "extra-focal image"
                    elif (channel=="B"):
                        wfsImgMap[sensorName].setImg(extraImg=wfsImg)

        return wfsImgMap

    def __searchDonutListId(self, donutList, starId):
        """
        
        Search the bright star ID in the donut list.
        
        Arguments:
            donutList {[list]} -- List of DonutImage object.
            starId {[int]} -- Star ID.
        
        Returns:
            [int] -- Index of donut image object with specific starId.
        """

        index = -1
        for ii in range(len(donutList)):
            if (donutList[ii].starId == int(starId)):
                index = ii
                break

        return index

    def getDonutMap(self, neighborStarMap, wfsImgMap, aFilter, doDeblending=False):
        """
        
        Get the donut map on each wavefront sensor (WFS).
        
        Arguments:
            neighborStarMap {[dict]} -- Information of neighboring stars and candidate stars with 
                                        the name of sensor as a dictionary.
            wfsImgMap {[dict]} --  Post-ISR image map.
            aFilter {[str]} -- Active filter type ("u", "g", "r", "i", "z", "y").
        
        Keyword Arguments:
            doDeblending {bool} -- Do the deblending or not. (default: {False})
        
        Returns:
            [dict] -- Donut image map.
        """

        # Get the corner wavefront sensor names
        wfsList = self.getWfsList()
        
        # Collect the donut images and put into the map/ dictionary
        donutMap = {}
        for sensorName in wfsImgMap.keys():

            # Get the abbraviated sensor name
            abbrevName = abbrevDectectorName(sensorName)

            # Configure the source processor
            self.sourProc.config(sensorName=abbrevName)

            # Get the bright star id list on specific sensor
            simobjIdList = list(neighborStarMap[sensorName].SimobjID.keys())

            # Get the defocal images: [intra, extra]
            defocalImgList = [wfsImgMap[sensorName].intraImg, wfsImgMap[sensorName].extraImg]

            for ii in range(len(simobjIdList)):
                
                # Get the single star map
                for jj in range(2):

                    ccdImg = defocalImgList[jj]

                    # Get the segment of image
                    if (ccdImg is not None):
                        singleSciNeiImg, allStarPosX, allStarPosY, magRatio, offsetX, offsetY = \
                                                        self.sourProc.getSingleTargetImage(ccdImg, 
                                                            neighborStarMap[sensorName], ii, aFilter)

                        # Get the single donut/ deblended image
                        imgDeblend = None
                        realcx = None
                        realcy = None
                        if ((len(magRatio) == 1 and doDeblending) or (not doDeblending)):
                            imgDeblend = singleSciNeiImg
                            realcx = allStarPosX[-1]
                            realcy = allStarPosY[-1]
                        # Do the deblending or not
                        elif (len(magRatio) == 2 and doDeblending):
                            imgDeblend, realcx, realcy = self.sourProc.doDeblending(singleSciNeiImg, 
                                                                  allStarPosX, allStarPosY, magRatio)

                        # Add the search algorithm here latter
                        if (len(magRatio) == 1):

                            # Calculate the center of mass first
                            realcy, realcx = center_of_mass(imgDeblend)

                        # Rotate the image if the sensor is the corner wavefront sensor
                        if sensorName in wfsList:

                            # Get the Euler angle
                            eulerZangle = round(self.sourProc.getEulerZinDeg(abbrevName))
                            
                            # Change the sign if the angle < 0
                            if (eulerZangle < 0):
                                eulerZangle += 360

                            # Do the rotation of matrix
                            numOfRot90 = eulerZangle//90
                            imgDeblend = np.flipud(np.rot90(np.flipud(imgDeblend), numOfRot90))

                        # Put the deblended image into the donut map
                        if (imgDeblend is not None):

                            if sensorName not in donutMap.keys():
                                donutMap[sensorName] = []

                            # Check the donut exists in the list or not
                            starId = simobjIdList[ii]
                            donutIndex = self.__searchDonutListId(donutMap[sensorName], starId)                             

                            # Create the donut object and put into the list if it is needed
                            if (donutIndex < 0):
                                donutImg = DonutImage(starId, realcx+offsetX, realcy+offsetY)
                                donutMap[sensorName].append(donutImg)

                                # Search for the donut index again
                                donutIndex = self.__searchDonutListId(donutMap[sensorName], starId)

                            # Get the donut image list
                            donutList = donutMap[sensorName]
                        
                            # Set the intra focal image
                            if (jj == 0):
                                donutList[donutIndex].setImg(intraImg=imgDeblend)
                            # Set the extra focal image
                            elif (jj == 1):
                                donutList[donutIndex].setImg(extraImg=imgDeblend)

        return donutMap

def plotDonutImg(donutMap, saveToDir=None, dpi=None):
    """
    
    Plot the donut image.
    
    Arguments:
        donutMap {[dict]} --  Donut image map.
    
    Keyword Arguments:
        saveToDir {[str]} -- Directory to save the images. (default: {None})
        dpi {[int]} -- The resolution in dots per inch. (default: {None})
    """

    intraType = "intra"
    extraType = "extra"

    for sensorName, donutList in donutMap.items():
        # Generate the image name
        imgTitle = abbrevDectectorName(sensorName) + "_DonutImg"

        # Collect all images and titles
        intraImgList = []
        extraImgList = []
        
        intraTitleList = []
        extraTitleList = []

        intraPixelXYList = []
        extraPixelXYList = []
        
        # Collect intra- and extra-focal donut images
        for donutImg in donutList:
            for ii in range(2):

                # Assign the image (0: intra, 1: extra)
                if (ii == 0):
                    img = donutImg.intraImg
                else:
                    img = donutImg.extraImg

                if (img is not None):

                    pixelXy = (donutImg.pixelX, donutImg.pixelY)

                    if (ii == 0):
                        intraImgList, intraTitleList, intraPixelXYList = _collectDonutImgList(
                                                intraImgList, intraTitleList, intraPixelXYList, 
                                                img, donutImg.starId, intraType, pixelXy)
                    else:
                        extraImgList, extraTitleList, extraPixelXYList = _collectDonutImgList(
                                                extraImgList, extraTitleList, extraPixelXYList, 
                                                img, donutImg.starId, extraType, pixelXy)

        # Decide the figure grid shape
        numOfRow = np.max([len(intraImgList), len(extraImgList)])

        if (len(intraImgList) == 0) or (len(extraImgList) == 0):
            numOfCol = 1
        else:
            numOfCol = 2

        gridShape = (numOfRow, numOfCol)

        # Plot the donut figure
        plt.figure()
        plt.suptitle(imgTitle)

        # Plot the intra-focal donut
        locOfCol = 0
        for ii in range(len(intraImgList)):
            _subPlot(plt, gridShape, (ii, locOfCol), intraImgList[ii], intraTitleList[ii], intraPixelXYList[ii])

        # Plot the extra-focal donut

        # Update the location of column if necessary
        if (numOfCol == 2):
            locOfCol = 1
        
        for ii in range(len(extraImgList)):
            _subPlot(plt, gridShape, (ii, locOfCol), extraImgList[ii], extraTitleList[ii], extraPixelXYList[ii])

        # Adjust the space between xlabel and title for neighboring sub-figures
        plt.tight_layout()

        # Save the file or not
        if (saveToDir is not None):
            # Generate the filepath
            imgType = ".png"
            imgFilePath = os.path.join(saveToDir, imgTitle+imgType)
            plt.savefig(imgFilePath, bbox_inches="tight", dpi=dpi)
            plt.close()
        else:
            plt.show()

def _subPlot(plt, gridShape, loc, img, aTitle, pixelXy):
    """
    
    Do the subplot of figure.
    
    Arguments:
        plt {[pyplot]} -- Plotting framework.
        gridShape {[tuple]} -- Shape of grid.
        loc {[tuple]} -- Location of subplot in grid.
        img {[ndarray]} -- Image of donut.
        aTitle {[str]} -- Title of subplot.
        pixelXy {[tuple]} -- Chip position of donut in (x, y).
    """

    # Chip position of donut
    pixelXy = np.round(pixelXy)
    pixelPos = "Pixel XY: (%d, %d)" % (pixelXy[0], pixelXy[1])

    # Decide the position of subplot
    ax = plt.subplot2grid(gridShape, loc)

    # Show the figure
    axPlot = ax.imshow(img, origin="lower")
    
    # Set the title
    ax.set_title(aTitle)

    # Set the x lavel
    ax.set_xlabel(pixelPos)

    # Set the colorbar
    plt.colorbar(axPlot, ax=ax)

def _collectDonutImgList(imgList, titleList, pixelXyList, img, starId, aType, pixelXy):
    """
    
    Collect the donut data in list.
    
    Arguments:
        imgList {[list]} -- List of image.
        titleList {[list]} -- List of title.
        pixelXyList {[list]} -- List of pixel XY.
        img {[ndarray]} -- Donut image.
        starId {[int]} -- Star Id.
        aType {[str]} -- Type of donut.
        pixelXy {[tuple]} -- Pixel position in (x, y).
    
    Returns:
        [list] -- List of image.
        [list] -- List of title.
        [list] -- List of pixel XY.
    """

    # Get the title
    aTitle = "_".join([str(starId), aType])    

    # Append the list
    imgList.append(img)
    titleList.append(aTitle)
    pixelXyList.append(pixelXy)

    return imgList, titleList, pixelXyList
    
if __name__ == "__main__":
    
    # Instintiate the components
    sourSelc = SourceSelector()
    dataCollector = WFDataCollector()
    isrWrapper = EimgIsrWrapper()
    sourProc = SourceProcessor()

    # Configurate the source selector
    # cameraType = "comcam"
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

    # Configurate the wfs data collector
    pathOfRawData = "../test/phosimOutput"
    destinationPath = "../test"
    butlerInputs = "../test"
    butlerOutputs = "../test"
    regisAdress = "../test/registry.sqlite3"
    dataCollector.config(pathOfRawData=pathOfRawData, destinationPath=destinationPath, 
                dbAdress=regisAdress, butlerInputs=butlerInputs, butlerOutputs=butlerOutputs)

    # Configurate the ISR wrapper
    isrWrapper.configWrapper(inputs=butlerInputs, outputs=butlerOutputs)

    # Configurate the source processor
    focalPlaneFolder = "../test"
    donutRadiusInPixel=63
    sourProc.config(donutRadiusInPixel=donutRadiusInPixel, folderPath2FocalPlane=focalPlaneFolder)

    # Initiate the WEP Controller
    wepCntlr = WEPController()
    wepCntlr.config(sourSelc=sourSelc, dataCollector=dataCollector, isrWrapper=isrWrapper, 
                    sourProc=sourProc)

    # Set the database address
    dbAdress = "../test/bsc.db3"

    # Do the query
    pointing = (0,0)
    cameraRotation = 0.0
    # skyInfoFilePath = "../test/phosimOutput/realComCam2/output/skyComCamInfo.txt"
    skyInfoFilePath = "../test/phosimOutput/realWfs/output/skyWfsInfo.txt"

    # camOrientation = "all"
    camOrientation = "corner"
    neighborStarMap, starMap, wavefrontSensors = wepCntlr.getTargetStarByFile(dbAdress, skyInfoFilePath, 
                                        pointing, cameraRotation, orientation=camOrientation, tableName="TempTable")

    # Import the PhoSim simulated image
    dataDirList = ["realComCam2/output/Extra", "realComCam2/output/Intra"]
    # for dataDir in dataDirList:
    #     wepCntlr.importPhoSimDataToButler(dataDir, atype="raw", overwrite=False)

    # Do the ISR
    extraObsId = 9005000
    intraObsId = 9005001
    obsIdList = [intraObsId, extraObsId]
    sensorNameList = list(starMap.keys())
    # for obsId in obsIdList:
    #     for sensorName in sensorNameList:
    #         wepCntlr.doISR(obsId, sensorName)

    # Get the wfs images
    # wfsImgMap = wepCntlr.getPostISRDefocalImgMap(sensorNameList, obsIdList=obsIdList)

    # Try the corner wavefront sensor
    # sensorNameList = wepCntlr.getWfsList()
    wfsDir = "realWfs/output"
    # cornerWfsImgMap = wepCntlr.getPostISRDefocalImgMap(sensorNameList, wfsDir=wfsDir)
    wfsImgMap = wepCntlr.getPostISRDefocalImgMap(sensorNameList, wfsDir=wfsDir)

    # Get the donut images
    donutMap = wepCntlr.getDonutMap(neighborStarMap, wfsImgMap, aFilter, doDeblending=False)

    # Check the donut
    for aKey, aItem in donutMap.items():
        donutList = aItem
        for ii in range(len(donutList)):
            print(aKey)
            print(donutList[ii].starId, donutList[ii].pixelX, donutList[ii].pixelY)

            if (donutList[ii].intraImg is not None):
                print(donutList[ii].intraImg.shape, np.sum(donutList[ii].intraImg))

            if (donutList[ii].extraImg is not None):
                print(donutList[ii].extraImg.shape, np.sum(donutList[ii].extraImg))

    # Plot the donut images
    saveToDir = "../test/donutImg"
    plotDonutImg(donutMap, saveToDir=saveToDir, dpi=None)
