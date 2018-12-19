import os, re, time, unittest
import numpy as np

import matplotlib
# Must be before importing matplotlib.pyplot or pylab!
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from astropy.io import fits
from scipy.ndimage.measurements import center_of_mass

from lsst.ts.wep.WFDataCollector import WFDataCollector
from lsst.ts.wep.SciIsrWrapper import SciIsrWrapper, getImageData
from lsst.ts.wep.SourceSelector import SourceSelector
from lsst.ts.wep.SourceProcessor import SourceProcessor, abbrevDectectorName
from lsst.ts.wep.WFEstimator import WFEstimator
from lsst.ts.wep.DefocalImage import DefocalImage, DonutImage
from lsst.ts.wep.Utility import getModulePath


class WEPController(object):

    def __init__(self):
        """
        
        Instantiate the WEP (wavefront estimation pipeline) controller.
        """
        
        self.sourSelc = None
        self.dataCollector = None
        self.isrWrapper = None
        self.sourProc = None
        self.wfsEsti = None

    def config(self, sourProc=None, dataCollector=None, isrWrapper=None, sourSelc=None, 
                wfsEsti=None):
        """
        
        Configurate the WEPController.
        
        Keyword Arguments:
            sourProc {[SourceProcessor]} -- Source Processor. (default: {None})
            dataCollector {[WFDataCollector]} -- Wavefront data collector. (default: {None})
            isrWrapper {[SciIsrWrapper]} -- ISR (instrument signature removal) wrapper. 
                                            (default: {None})
            sourSelc {[SourceSelector]} -- Source selector. (default: {None})
            wfsEsti {[WFEstimator]} -- Wavefront estimator. (default: {None})
        """

        self.__setVar(sourSelc, "sourSelc")        
        self.__setVar(dataCollector, "dataCollector")
        self.__setVar(isrWrapper, "isrWrapper")
        self.__setVar(sourProc, "sourProc")
        self.__setVar(wfsEsti, "wfsEsti")

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

    def ingestSimImages(self, fitsFileArg=None, dataDir=None, atype="raw", overwrite=False):
        """
        
        Import the PhoSim simulated data to match with the data butler to use. This means the 
        registry.sqlite3 repo will be inserted with the meta data if necessary.
        
        Keyword Arguments:
            fitsFileArg {[str]} -- Fits file argument. This is for DM cmd task. (default: {None})
            dataDir {[str]} -- PhoSim FITS data directory. (default: {None})
            atype {[str]} -- Dataset type. (default: {"raw"})
            overwrite {[boolean]} -- Overwrite the existed files or not. (default: {False})
        
        Raises:
            ValueError -- Not allowed type ("raw", "bias", "dark", "flat").
        """

        if (fitsFileArg is not None):
            self.dataCollector.ingestSimImages(fitsFileArg=fitsFileArg)

        else:

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
            self.dataCollector.ingestSimImages(dataDir, obsId=obsId, aFilter=aFilter, atype=atype, 
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
            if (isinstance(self.isrWrapper, SciIsrWrapper)):
                self.isrWrapper.doISR(visit, snap, raft, sensor)
            else:
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

    def getPostISRDefocalImgMap(self, sensorNameList, obsIdList=None, wfsDir=None, snap=0, expInDmCoor=False):
        """
        
        Get the post-ISR defocal image map.
        
        Arguments:
            sensorNameList {[list]} -- List of sensor name which is in the canonical form.
        
        Keyword Arguments:
            obsIdList {[list]} -- Observation Id list in [intraObsId, extraObsId]. (default: {None})
            wfsDir {[str]} -- Directory to wavefront sensor image data. (default: {None})
            snap {[int]} -- Snap number (default: {0})
            expInDmCoor {[bool]} -- Exposure image is in DM coordinate system. If True, this function will 
                                    rotate the exposure image to camera coordinate. This only works for 
                                    LSST FAM at this moment.
        
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

                    # Get the exposure image in ndarray
                    img = getImageData(img)

                    # Change the image to camera coordinate
                    if (expInDmCoor):
                        img = np.rot90(img.copy(), k=3)

                    imgList.append(img)

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

    def getDonutMap(self, neighborStarMap, wfsImgMap, aFilter, doDeblending=False, sglDonutOnly=False):
        """
        
        Get the donut map on each wavefront sensor (WFS).
        
        Arguments:
            neighborStarMap {[dict]} -- Information of neighboring stars and candidate stars with 
                                        the name of sensor as a dictionary.
            wfsImgMap {[dict]} --  Post-ISR image map.
            aFilter {[str]} -- Active filter type ("u", "g", "r", "i", "z", "y").
        
        Keyword Arguments:
            doDeblending {[bool]} -- Do the deblending or not. (default: {False})
            sglDonutOnly{[bool]} -- Only consider the single donut based on the bright star catalog. 
                                    (default: {False})
        
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

                        # Check the single donut or not based on the bright star catalog only
                        # This method should be updated in the future
                        if (sglDonutOnly):
                            if (len(magRatio) != 1):
                                continue

                        # Get the single donut/ deblended image
                        if (len(magRatio) == 1) or (not doDeblending):
                            imgDeblend = singleSciNeiImg

                            if (len(magRatio) == 1):
                                realcx, realcy = searchDonutPos(imgDeblend)
                            else:                               
                                realcx = allStarPosX[-1]
                                realcy = allStarPosY[-1]
                        # Do the deblending or not
                        elif (len(magRatio) == 2 and doDeblending):
                            imgDeblend, realcx, realcy = self.sourProc.doDeblending(singleSciNeiImg, 
                                                                  allStarPosX, allStarPosY, magRatio)
                            # Update the magnitude ratio
                            magRatio = [1]

                        else:
                            continue

                        # Extract the image
                        if (len(magRatio) == 1):
                            x0 = np.floor(realcx-self.wfsEsti.sizeInPix/2).astype("int")
                            y0 = np.floor(realcy-self.wfsEsti.sizeInPix/2).astype("int")
                            imgDeblend = imgDeblend[y0:y0+self.wfsEsti.sizeInPix, 
                                                    x0:x0+self.wfsEsti.sizeInPix]

                        # Rotate the image if the sensor is the corner wavefront sensor
                        if sensorName in wfsList:

                            # Get the Euler angle
                            eulerZangle = round(self.sourProc.getEulerZinDeg(abbrevName))
                            
                            # Change the sign if the angle < 0
                            while (eulerZangle < 0):
                                eulerZangle += 360

                            # Do the rotation of matrix
                            numOfRot90 = eulerZangle//90
                            imgDeblend = np.flipud(np.rot90(np.flipud(imgDeblend), numOfRot90))

                        # Put the deblended image into the donut map
                        if sensorName not in donutMap.keys():
                            donutMap[sensorName] = []

                        # Check the donut exists in the list or not
                        starId = simobjIdList[ii]
                        donutIndex = self.__searchDonutListId(donutMap[sensorName], starId)                             

                        # Create the donut object and put into the list if it is needed
                        if (donutIndex < 0):

                            # Calculate the field X, Y
                            pixelX = realcx+offsetX
                            pixelY = realcy+offsetY
                            fieldX, fieldY = self.sourProc.camXYtoFieldXY(pixelX, pixelY)

                            # Instantiate the DonutImage class
                            donutImg = DonutImage(starId, pixelX, pixelY, fieldX, fieldY)
                            donutMap[sensorName].append(donutImg)

                            # Search for the donut index again
                            donutIndex = self.__searchDonutListId(donutMap[sensorName], starId)

                        # Get the donut image list
                        donutList = donutMap[sensorName]

                        # Take the absolute value for images, which might contain the 
                        # negative value after the ISR correction. This happens for the 
                        # amplifier images.
                        imgDeblend = np.abs(imgDeblend)
                    
                        # Set the intra focal image
                        if (jj == 0):
                            donutList[donutIndex].setImg(intraImg=imgDeblend)
                        # Set the extra focal image
                        elif (jj == 1):
                            donutList[donutIndex].setImg(extraImg=imgDeblend)

        return donutMap

    def genMasterImgSglCcd(self, sensorName, donutImgList, zcCol=np.zeros(22)):
        """
        
        Generate the master donut image on signle CCD.
        
        Arguments:
            sensorName {[str]} -- Sensor name.
            donutImgList {[list]} -- List of donut images.
        
        Keyword Arguments:
            zcCol {[ndarray]} -- Coefficients of wavefront (z1-z22). (default: {np.zeros(22)})
        
        Returns:
            [DefocalImage] -- Master donut image.
        """

        # Configure the source processor
        abbrevName = abbrevDectectorName(sensorName)
        self.sourProc.config(sensorName=abbrevName)

        intraProjImgList = []
        extraProjImgList = []

        for donutImg in donutImgList:

            # Get the field x, y
            pixelX = donutImg.pixelX
            pixelY = donutImg.pixelY

            fieldX, fieldY = self.sourProc.camXYtoFieldXY(pixelX, pixelY)
            fieldXY = (fieldX, fieldY)

            # Set the image
            if (donutImg.intraImg is not None):

                # Get the projected image
                projImg = self.__getProjImg(fieldXY, donutImg.intraImg, 
                                            self.wfsEsti.ImgIntra.INTRA, zcCol)

                # Collect the projected donut
                intraProjImgList.append(projImg)

            if (donutImg.extraImg is not None):
                
                # Get the projected image
                projImg = self.__getProjImg(fieldXY, donutImg.extraImg, 
                                            self.wfsEsti.ImgExtra.EXTRA, zcCol)

                # Collect the projected donut
                extraProjImgList.append(projImg)

        # Generate the master donut
        stackIntraImg = self.__stackImg(intraProjImgList)
        stackExtraImg = self.__stackImg(extraProjImgList)

        # Put the master donut to donut map
        masterDonut = DonutImage(0, None, None, 0, 0, intraImg=stackIntraImg, 
                                    extraImg=stackExtraImg)

        return masterDonut

    def generateMasterImg(self, donutMap, zcCol=np.zeros(22)):
        """
        
        Generate the master donut image map.
        
        Arguments:
            donutMap {[dict]} -- Donut image map.
        
        Keyword Arguments:
            defocalDisInMm {float} -- Defocal distance in mm. (default: {1})
            zcCol {[ndarray]} -- Coefficients of wavefront (z1-z22) in m. 
                                 (default: {np.zeros(22)})
        
        Returns:
            [dict] -- Master donut image map.
        """

        masterDonutMap = {}
        for sensorName, donutImgList in donutMap.items():

            # Get the master donut on single CCD
            masterDonut = self.genMasterImgSglCcd(sensorName, donutImgList, zcCol=zcCol)

            # Put the master donut to donut map
            masterDonutMap[sensorName] = [masterDonut]

        return masterDonutMap

    def __stackImg(self, imgList):
        """
        
        Stack the images to generate the master image.
        
        Arguments:
            imgList {[list]} -- Image list.
        
        Returns:
            [ndarray] -- Stacked image.
        """

        if (len(imgList) == 0):
            stackImg = None
        else:
            # Get the minimun image dimension
            dimXlist = []
            dimYlist = []
            for img in imgList:
                dy, dx = img.shape
                dimXlist.append(dx)
                dimYlist.append(dy)

            dimX = np.min(dimXlist)
            dimY = np.min(dimYlist)

            deltaX = dimX//2
            deltaY = dimY//2

            # Stack the image by summation directly
            stackImg = np.zeros([dimY, dimX])
            for img in imgList:
                dy, dx = img.shape
                cy = int(dy/2)
                cx = int(dx/2)
                stackImg += img[cy-deltaY:cy+deltaY, cx-deltaX:cx+deltaX]

        return stackImg

    def __getProjImg(self, fieldXY, defocalImg, aType, zcCol):
        """
        
        Get the projected image on the pupil.
        
        Arguments:
            fieldXY {[tuple]} -- Position of donut on the focal plane in degree for intra- and 
                                 extra-focal images.
            defocalImg {[ndarray]} -- Defocal image.
            aType {[str]} -- Defocal type.
            zcCol {[ndarray]} -- Coefficients of wavefront (z1-z22).
        
        Returns:
            [ndarray] -- Projected image.
        """
        
        # Set the image
        self.wfsEsti.setImg(fieldXY, image=defocalImg, defocalType=aType)

        # Get the distortion correction (offaxis)
        offAxisCorrOrder = self.wfsEsti.algo.parameter["offAxisPolyOrder"]
        instDir = os.path.dirname(self.wfsEsti.inst.filename)
        if (aType == self.wfsEsti.ImgIntra.INTRA):
            img = self.wfsEsti.ImgIntra
        elif (aType == self.wfsEsti.ImgExtra.EXTRA):
            img = self.wfsEsti.ImgExtra
        img.getOffAxisCorr(instDir, offAxisCorrOrder)

        # Do the image cocenter
        img.imageCoCenter(self.wfsEsti.inst)

        # Do the compensation/ projection
        img.compensate(self.wfsEsti.inst, self.wfsEsti.algo, zcCol, self.wfsEsti.opticalModel)

        # Return the projected image
        return img.image

    def calcWfErr(self, donutMap):
        """
        
        Calculate the wavefront error in annular Zernike polynomials (z4-z22).
        
        Arguments:
            donutMap {[dict]} -- Donut image map.
        
        Returns:
            [dict] -- Donut image map with calculated wavefront error.
        """

        for sensorName, donutList in donutMap.items():

            for ii in range(len(donutList)):

                # Get the intra- and extra-focal donut images

                # Check the sensor is the corner WFS or not. Only consider "A"
                # Intra: C0 -> A; Extra: C1 -> B

                # Look for the intra-focal image
                if sensorName.endswith("A"):
                    intraDonut = donutList[ii]

                    # Get the extra-focal sensor name
                    extraFocalSensorName = sensorName.replace("A", "B")

                    # Get the donut list of extra-focal sensor
                    extraDonutList = donutMap[extraFocalSensorName]
                    if (ii < len(extraDonutList)):
                        extraDonut = extraDonutList[ii]
                    else:
                        continue

                # Pass the extra-focal image
                elif sensorName.endswith("B"):
                    continue
                # Scientific sensor
                else:
                    intraDonut = extraDonut = donutList[ii]

                # Calculate the wavefront error

                # Get the field X, Y position
                intraFieldXY = (intraDonut.fieldX, intraDonut.fieldY)
                extraFieldXY = (extraDonut.fieldX, extraDonut.fieldY)

                # Get the defocal images
                intraImg = intraDonut.intraImg
                extraImg = extraDonut.extraImg

                # Calculate the wavefront error
                zer4UpNm = self.calcSglWfErr(intraImg, extraImg, intraFieldXY, extraFieldXY)

                # Put the value to the donut image
                intraDonut.setWfErr(zer4UpNm)
                extraDonut.setWfErr(zer4UpNm)

        return donutMap

    def calcSglWfErr(self, intraImg, extraImg, intraFieldXY, extraFieldXY):
        """
        
        Calculate the wavefront error in annular Zernike polynomials (z4-z22) for 
        single donut.
        
        Arguments:
            intraImg {[ndarray]} -- Intra-focal donut image.
            extraImg {[ndarray]} -- Extra-focal donut image.
            intraFieldXY {[tuple]} -- Field x, y in degree of intra-focal donut image.
            extraFieldXY {[tuple]} -- Field x, y in degree of extra-focal donut image.
        
        Returns:
            [ndarray] -- Coefficients of Zernike polynomials (z4 - z22) in nm.
        """

        # Set the images
        self.wfsEsti.setImg(intraFieldXY, image=intraImg, 
                            defocalType=self.wfsEsti.ImgIntra.INTRA)
        self.wfsEsti.setImg(extraFieldXY, image=extraImg, 
                            defocalType=self.wfsEsti.ImgExtra.EXTRA)

        # Reset the wavefront estimator
        self.wfsEsti.reset()

        # Calculate the wavefront error
        zer4UpNm = self.wfsEsti.calWfsErr()

        return zer4UpNm

    def calcSglAvgWfErr(self, donutImgList):
        """
        
        Calculate the average of wavefront error on single CCD.
        
        Arguments:
            donutImgList {[list]} -- List of donut images.
        
        Returns:
            [ndarray] -- Average of wavefront error in nm.
        """

        # Calculate the weighting of donut image
        weightingRatio = calcWeiRatio(donutImgList)

        # Calculate the mean wavefront error (z4 - z22)
        numOfZk = self.wfsEsti.algo.parameter["numTerms"] - 3
        avgErr = np.zeros(numOfZk)
        for ii in range(len(donutImgList)):
            donutImg = donutImgList[ii]

            # Get the zer4UpNm
            if (donutImg.zer4UpNm is not None):
                zer4UpNm = donutImg.zer4UpNm
            else:
                zer4UpNm = 0

            avgErr = avgErr + weightingRatio[ii]*zer4UpNm

        return avgErr

    def sortDonut(self, donutList):

        # This function is to sort the donut images in list from high S/N to low.
        pass

def calcWeiRatio(donutImgList):
    """
    
    Calculate the weighting ratio of donut image in the list.
    
    Arguments:
        donutImgList {[list]} -- List of donut images.
    
    Returns:
        [ndarray] -- Array of Weighting ratio of image donuts.
    """

    # Weighting of donut image. Use the simple average at this moment.
    # Need to consider the S/N and other factors in the future

    # Check the available zk and give the ratio
    weightingRatio = []
    for donutImg in donutImgList:
        if (donutImg.zer4UpNm is not None):
            weightingRatio.append(1)
        else:
            weightingRatio.append(0)

    # Do the normalization
    weightingRatio = np.array(weightingRatio)
    weightingRatio = weightingRatio/np.sum(weightingRatio)

    return weightingRatio

def searchDonutPos(img):
    """
    
    Search the position of donut on image.
    
    Arguments:
        img {[ndarray]} -- Donut image.
    
    Returns:
        [float] -- x position of donut center in pixel.
        [float] -- y position of donut center in pixel.
    """

    # Search the donut position by the center of mass
    # Need to update this method to the more robust one such as the convolution
    realcy, realcx = center_of_mass(img)

    return realcx, realcy

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
    
class WEPControllerTest(unittest.TestCase):

    """ 
    Test the function of WEPController.
    """

    def setUp(self):

        # Get the path of module
        self.modulePath = getModulePath()

        # Instantiate the WEP controller
        self.wepCntlr = WEPController()

    def testCornerWfsFunction(self):

        # Test to get the list of corner wavefront sensors
        wfsList = self.wepCntlr.getWfsList()
        self.assertEqual(len(wfsList), 8)

        # Instintiate the components
        sourSelc = SourceSelector()
        dataCollector = WFDataCollector()
        sourProc = SourceProcessor()

        instruFolderPath = os.path.join(self.modulePath, "algoData", "cwfs", "instruData")
        algoFolderPath = os.path.join(self.modulePath, "algoData", "cwfs", "algo")
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
        pathOfRawData = os.path.join(self.modulePath, "test", "phosimOutput")
        destinationPath = butlerInputs = butlerOutputs = os.path.join(self.modulePath, "test")
        dataCollector.config(pathOfRawData=pathOfRawData, destinationPath=destinationPath)

        # Configurate the source processor
        focalPlaneFolder = os.path.join(self.modulePath, "test")
        sourProc.config(donutRadiusInPixel=starRadiusInPixel, folderPath2FocalPlane=focalPlaneFolder, 
                        pixel2Arcsec=0.2)

        # Configurate the wavefront estimator
        defocalDisInMm = None
        
        # Size of donut in pixel for corner WFS
        sizeInPix = 120
        wfsEsti.config(solver="exp", instName=cameraType, opticalModel="offAxis", 
                        defocalDisInMm=defocalDisInMm, sizeInPix=sizeInPix)

        # Configurate the WEP controller
        self.wepCntlr.config(sourSelc=sourSelc, dataCollector=dataCollector, 
                            sourProc=sourProc, wfsEsti=wfsEsti)

        # Test the configuration
        self.assertTrue(isinstance(self.wepCntlr.wfsEsti, WFEstimator))

        # Get the target stars by file

        # Set the database address
        dbAdress = os.path.join(self.modulePath, "test", "bsc.db3")

        # Do the query
        pointing = (0,0)
        cameraRotation = 0.0
        skyInfoFilePath = os.path.join(self.modulePath, "test", "phosimOutput", "realWfs", "output", 
                                       "skyWfsInfo.txt")

        camOrientation = "corner"
        neighborStarMap, starMap, wavefrontSensors = self.wepCntlr.getTargetStarByFile(dbAdress, 
                                                        skyInfoFilePath, pointing, cameraRotation, 
                                                orientation=camOrientation, tableName="TempTable")
        self.assertEqual(len(starMap), 8)

        starData = starMap["R:0,0 S:2,2,A"]
        self.assertEqual(len(starData.SimobjID), 2)

        # Get the available sensor name list
        sensorNameList = list(starMap.keys())

        # Get the eimage
        wfsDir = os.path.join("realWfs", "output")
        wfsImgMap = self.wepCntlr.getPostISRDefocalImgMap(sensorNameList, wfsDir=wfsDir)

        wfsImg = wfsImgMap["R:0,0 S:2,2,A"]
        self.assertNotEqual(np.sum(wfsImg.intraImg), None)
        self.assertEqual(wfsImg.extraImg, None)

        # Get the donut map
        donutMap = self.wepCntlr.getDonutMap(neighborStarMap, wfsImgMap, aFilter, 
                                            doDeblending=False, sglDonutOnly=True)
        
        donutList = donutMap["R:0,0 S:2,2,A"]
        self.assertEqual(len(donutList), 2)

        donutImg = donutList[0]
        self.assertNotEqual(np.sum(donutImg.intraImg), None)
        self.assertEqual(donutImg.extraImg, None)
        self.assertEqual(donutImg.starId, 6)
        self.assertEqual(int(donutImg.pixelX), 506)
        self.assertEqual(int(donutImg.pixelY), 1008)

        # Calculate the wavefront error for the individual donut
        partDonutMap = dict()
        partDonutMap["R:0,0 S:2,2,A"] = donutMap["R:0,0 S:2,2,A"]
        partDonutMap["R:0,0 S:2,2,B"] = donutMap["R:0,0 S:2,2,B"]
        
        partDonutMap = self.wepCntlr.calcWfErr(partDonutMap)
        
        donutList = partDonutMap["R:0,0 S:2,2,A"]
        donutImg = donutList[0]
        self.assertEqual(len(donutImg.zer4UpNm), 19)

        # Test the weighting ratio
        weightingRatio = calcWeiRatio(donutList)
        self.assertEqual(np.sum(weightingRatio), 1)
        self.assertEqual(weightingRatio[0], 0.5)

        # Test to calculate the average wavefront error
        avgErr = self.wepCntlr.calcSglAvgWfErr(donutList)
        ans = donutList[0].zer4UpNm*weightingRatio[0] + donutList[1].zer4UpNm*weightingRatio[1]
        self.assertEqual(np.sum(avgErr), np.sum(ans))


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
