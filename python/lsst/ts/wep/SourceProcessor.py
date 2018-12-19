import os, re
import numpy as np

from lsst.ts.wep.deblend.BlendedImageDecorator import BlendedImageDecorator
from lsst.ts.wep.isr.changePhoSimInstrument import readData
from lsst.ts.wep.SourceSelector import SourceSelector
from lsst.ts.wep.SciIsrWrapper import poltExposureImage


class SourceProcessor(object):

    def __init__(self):
        """
        
        Initialize the SourceProcessor class.
        """

        self.sensorName = None
        self.donutRadiusInPixel = None

        self.sensorFocaPlaneInDeg = None
        self.sensorFocaPlaneInUm = None
        self.sensorDimList = None
        self.sensorEulerRot = None

        self.pixel2Arcsec = None
        self.pixel2um = None

        self.blendedImageDecorator = BlendedImageDecorator()

    def config(self, sensorName=None, donutRadiusInPixel=None, folderPath2FocalPlane=None, 
                pixel2Arcsec=0.2, pixel2um=10.0):
        """
        
        Do the configuration.
        
        Keyword Arguments:
            sensorName {[str]} -- Sensor name. (default: {None})
            donutRadiusInPixel {[float]} -- Donut radius in pixel. (default: {None})
            folderPath2FocalPlane {[str]} -- Path to the directory of focal plane data 
                                             ("focalplanelayout.txt"). (default: {None})
            pixel2Arcsec {[float]} -- [Pixel to arcsec. (default: {0.2})
            pixel2um {[float]} -- Pixel to um. (default: {10.0})
        """

        # Give the sensor name
        if (sensorName is not None):
            self.sensorName = sensorName

        # Give the donut radius in pixel
        if (donutRadiusInPixel is not None):
            self.donutRadiusInPixel = donutRadiusInPixel

        # Read the focal plane data
        if (folderPath2FocalPlane is not None):
            self.__readFocalPlane(folderPath2FocalPlane, pixel2Arcsec=pixel2Arcsec)

        # Set the unit of pixel to arcsec
        self.pixel2Arcsec = pixel2Arcsec

        # Set the unit of pixel to micron
        self.pixel2um = pixel2um

    def __readFocalPlane(self, folderPath, pixel2Arcsec, fileName="focalplanelayout.txt"):
        """

        Read the focal plane data used in PhoSim to get the ccd dimension and fieldXY in chip 
        center.

        Arguments:
            folderPath {[str]} -- Directory of focal plane file.
            pixel2Arcsec {float} -- Pixel to arcsec.

        Keyword Arguments:
            fileName {[str]} -- Filename of focal plane. (default: {"focalplanelayout.txt"})
        """

        # Read the focal plane data by the delegation
        ccdData = readData(folderPath, fileName, "fieldCenter")

        # Collect the focal plane data
        sensorFocaPlaneInDeg = {}
        sensorFocaPlaneInUm = {}
        sensorDimList = {}
        for akey, aitem in ccdData.items():

            # Consider the x-translation in corner wavefront sensors
            aitem = self.__shiftCenterWfs(akey, aitem)

            # Change the unit from um to degree
            fieldX = float(aitem[0])/float(aitem[2])*pixel2Arcsec/3600
            fieldY = float(aitem[1])/float(aitem[2])*pixel2Arcsec/3600

            # Get the data
            sensorFocaPlaneInDeg.update({akey: (fieldX, fieldY)})
            sensorFocaPlaneInUm.update({akey: (float(aitem[0]), float(aitem[1]))})
            sensorDimList.update({akey: (int(aitem[3]), int(aitem[4]))})

        # Assign the values
        self.sensorDimList = sensorDimList
        self.sensorFocaPlaneInDeg = sensorFocaPlaneInDeg
        self.sensorFocaPlaneInUm = sensorFocaPlaneInUm
        self.sensorEulerRot = readData(folderPath, fileName, "eulerRot")

    def __shiftCenterWfs(self, sensorName, focalPlaneData):
        """

        Get the fieldXY of center of wavefront sensors. The input data is the center of 
        combined chips (C0+C1). The layout is shown in the following:

        R04_S20              R44_S00
        --------           -----------       /\ +y
        |  C0  |           |    |    |        |
        |------|           | C1 | C0 |        |
        |  C1  |           |    |    |        |
        --------           -----------        -----> +x

        R00_S22              R40_S02
        -----------          --------
        |    |    |          |  C1  |
        | C0 | C1 |          |------|
        |    |    |          |  C0  |
        -----------          --------

        Arguments:
            sensorName {[str]} -- Sensor name.
            focalPlaneData {[list]} -- Data of focal plane: x position (microns), y position 
                                       (microns), pixel size (microns), number of x pixels, 
                                       number of y pixels.
        Returns:
            [list] -- Updated focal plane data.
        """

        # Consider the x-translation in corner wavefront sensors
        tempX = None
        tempY = None

        if sensorName in ("R44_S00_C0", "R00_S22_C1"):
            # Shift center to +x direction
            tempX = float(focalPlaneData[0]) + float(focalPlaneData[3])/2*float(focalPlaneData[2])
        elif sensorName in ("R44_S00_C1", "R00_S22_C0"):
            # Shift center to -x direction
            tempX = float(focalPlaneData[0]) - float(focalPlaneData[3])/2*float(focalPlaneData[2])
        elif sensorName in ("R04_S20_C1", "R40_S02_C0"):
            # Shift center to -y direction
            tempY = float(focalPlaneData[1]) - float(focalPlaneData[3])/2*float(focalPlaneData[2])
        elif sensorName in ("R04_S20_C0", "R40_S02_C1"):
            # Shift center to +y direction
            tempY = float(focalPlaneData[1]) + float(focalPlaneData[3])/2*float(focalPlaneData[2])

        # Replace the value by the shifted one
        if (tempX is not None):
            focalPlaneData[0] = str(tempX)
        elif (tempY is not None):
            focalPlaneData[1] = str(tempY)

        # Return the center position of wave front sensor
        return focalPlaneData

    def camXYtoFieldXY(self, pixelX, pixelY):
        """

        Get the field X, Y of the pixel postion in CCD. It is noted that the wavefront 
        sensors will do the counter-clockwise rotation as the following based on the euler 
        angle:

        R04_S20              R44_S00
        O-------           -----O----O       /\ +y
        |  C0  |           |    |    |        |
        O------|           | C1 | C0 |        |
        |  C1  |           |    |    |        |
        --------           -----------        O----> +x

        R00_S22              R40_S02
        -----------          --------
        |    |    |          |  C1  |
        | C0 | C1 |          |------O
        |    |    |          |  C0  |
        O----O-----          -------O

        Arguments:
            pixelX {[float]} -- Pixel x on camera coordinate.
            pixelY {[float]} -- Pixel y on camera coordinate.

        Returns:
            [float] -- Field X, Y in degree.
        """

        # Get the field X, Y of sensor's center
        fieldXc, fieldYc = self.sensorFocaPlaneInDeg[self.sensorName]

        # Get the center pixel position
        pixelXc, pixelYc = self.sensorDimList[self.sensorName]
        pixelXc = pixelXc/2
        pixelYc = pixelYc/2

        # Calculate the delta x and y in degree
        deltaX = (pixelX-pixelXc)*self.pixel2Arcsec/3600.0
        deltaY = (pixelY-pixelYc)*self.pixel2Arcsec/3600.0

        # Calculate the transformed coordinate in degree.
        fieldX, fieldY = self.__rotCam2FocalPlane(self.sensorName, fieldXc, fieldYc, 
                                                  deltaX, deltaY)

        return fieldX, fieldY

    def focalPlaneXY2CamXY(self, xInUm, yInUm):
        """
        
        Get the x, y position on camera plane from the focal plane position.
        
        Arguments:
            xInUm {[float]} -- Position x on focal plane in um.
            yInUm {[float]} -- Position y on focal plane in um.

        Returns:
            [float] -- Pixel x, y position on camera plane.
        """

        # Get the central position of sensor in um
        xc, yc = self.sensorFocaPlaneInUm[self.sensorName]

        # Get the center pixel position
        pixelXc, pixelYc = self.sensorDimList[self.sensorName]
        pixelXc = pixelXc/2
        pixelYc = pixelYc/2

        # Calculate the delta x and y in pixel
        deltaX = (xInUm-xc)/self.pixel2um
        deltaY = (yInUm-yc)/self.pixel2um

        # Calculate the transformed coordinate
        pixelX, pixelY = self.__rotCam2FocalPlane(self.sensorName, pixelXc, pixelYc, deltaX, 
                                                    deltaY, counterClockWise=False)

        return pixelX, pixelY

    def __rotCam2FocalPlane(self, sensorName, centerX, centerY, deltaX, deltaY, 
                            counterClockWise=True):
        """
        
        Do the rotation from camera coordinate to focal plane coordinate or vice versa.
        
        Arguments:
            sensorName {[str]} -- Sensor name.
            centerX {[float]} -- CCD center X.
            centerY {[float]} -- CCD center Y.
            deltaX {[float]} -- Delta X from the CCD's center.
            deltaY {[float]} -- Delta Y from the CCD's center.
        
        Keyword Arguments:
            counterClockWise {bool} -- Direction of rotation: counter-clockwise or 
                                        clockwise. (default: {True})
        
        Returns:
            [float] -- Transformed X, Y position.
        """

        # Get the euler angle in z direction (only consider the z rotatioin at this moment)
        eulerZ = round(self.getEulerZinDeg(sensorName))

        # Change the unit to radian
        eulerZ = eulerZ/180.0*np.pi

        # Counter-clockwise or clockwise rotation
        if (not counterClockWise):
            eulerZ = -eulerZ

        # Calculate the new x, y by the rotation. This is important for wavefront sensor.
        newX = centerX + np.cos(eulerZ)*deltaX - np.sin(eulerZ)*deltaY
        newY = centerY + np.sin(eulerZ)*deltaX + np.cos(eulerZ)*deltaY

        return newX, newY

    def getEulerZinDeg(self, abbrevName):
        """
        
        Get the Euler Z angle of sensor in degree.
        
        Arguments:
            abbrevName {[str]} -- Abbreviated sensor name.
        
        Returns:
            [float] -- Euler Z angle in degree.
        """

        return float(self.sensorEulerRot[abbrevName][0])

    def dmXY2CamXY(self, pixelDmX, pixelDmY):
        """

        Transform the pixel x, y from DM library to camera to use. Camera coordinate is defined
        in LCA-13381. Define camera coordinate (x, y) and DM coordinate (x', y'), then the relation
        is dx' = -dy, dy' = dx.

         O---->y
         |
         |   ----------------------
         \/ |                      |   (x', y') = (200, 500) => (x, y) = (-500, 200) -> (3500, 200)
         x  |                      |
            |4000                  |
        y'  |                      |
         /\ |       4072           |
         |  |----------------------
         |
         O-----> x'

        Arguments:
            pixelDmX {[float]} -- Pixel x defined in DM coordinate.
            pixelDmY {[float]} -- Pixel y defined in DM coordinate.

        Returns:
            [float] -- Pixel x, y defined in camera coordinate based on LCA-13381.
        """

        # Get the CCD dimension
        dimX, dimY = self.sensorDimList[self.sensorName]

        # Calculate the transformed coordinate
        pixelCamX = dimX-pixelDmY
        pixelCamY = pixelDmX

        return pixelCamX, pixelCamY

    def camXY2DmXY(self, pixelCamX, pixelCamY):
        """
        
        Transform the pixel x, y from camera coordinate to DM coordinate. Camera coordinate is 
        defined in LCA-13381.
        
        Arguments:
            pixelCamX {[float]} -- Pixel x defined in Camera coordinate based on LCA-13381.
            pixelCamY {[float]} -- Pixel y defined in Camera coordinate based on LCA-13381.
        
        Returns:
            [float] -- Pixel x, y defined in DM coordinate.
        """
        
        # Get the CCD dimension
        dimX, dimY = self.sensorDimList[self.sensorName]

        # Calculate the transformed coordinate
        pixelDmX = pixelCamY
        pixelDmY = dimX-pixelCamX

        return pixelDmX, pixelDmY

    def evalVignette(self, fieldX, fieldY, distanceToVignette=1.75):
        """
        
        Evaluate the donut is vignetted or not by comparing the donut's distance to center with
        a reference value.
        
        Arguments:
            fieldX {[float]} -- Field x in degree.
            fieldY {[float]} -- Field y in degree.
        
        Keyword Arguments:
            distanceToVignette {float} -- Reference to be the vignetting. Use the half of field 
                                            of view as a initial guess. (default: {1.75})
        
        Returns:
            [type] -- [description]
        """
        
        # The donut is vignetted or not.
        isVignette = False

        # Calculate the distance to center in degree to judge the donut is vignetted or not.
        fldr = np.sqrt(fieldX**2 + fieldY**2)
        if (fldr >= distanceToVignette):
            isVignette = True 

        return isVignette

    def doDeblending(self, blendedImg, allStarPosX, allStarPosY, magRatio):
        """
        
        Do the deblending. It is noted that the algorithm now is only for one bright star and 
        one neighboring star.
        
        Arguments:
            blendedImg {[float]} -- Blended image.
            allStarPosX {[float]} -- Star's position x in pixel. The final one is the bright star.
            allStarPosY {[float]} -- Star's position y in pixel. The final one is the bright star.
            magRatio {[float]} -- Star's magnitude compared with the bright star.
        
        Returns:
            [float] -- Deblended image.
            [float] -- Pixel x, y of bright star.
        
        Raises:
            ValueError -- The inputs are not one bright star + one neighboring star.
        """

        # Check there is only one bright star and one neighboring star. This is the limit of 
        # deblending algorithm now.
        if (len(magRatio) != 2):
            raise ValueError("Deblending can only handle one bright star and one neighboring Star now.")

        # Set the image for the deblending
        self.blendedImageDecorator.setImg(image=blendedImg)

        # Do the deblending
        imgDeblend, realcx, realcy = self.blendedImageDecorator.deblendDonut([allStarPosX[0], 
                                                                 allStarPosY[0]], magRatio[0])

        return imgDeblend, realcx, realcy

    def getSingleTargetImage(self, ccdImg, neighboringStarMapOnSingleSensor, index, aFilter):
        """

        Get the image of single scientific target and related neighboring stars.

        Arguments:
            ccdImg {[float]} -- CCD image.
            neighboringStarMapOnSingleSensor {[dict]} -- Neighboring star map.
            index {[int]} -- Index of science target star in neighboring star map.
            aFilter {[str]} -- Active filter type

        Returns:
            [float] -- Ccd image of target stars.
            [float] -- Star positions in x, y.
            [float] -- Star magnitude ratio compared with the bright star.
            [float] -- Offset x, y from the origin of target star image to the origin of 
                        CCD image.

        Raises:
            ValueError -- Science star index is out of the neighboring star map.
        """
    
        # Get the target star position
        if (index >= len(neighboringStarMapOnSingleSensor.SimobjID)):
            raise ValueError("Index is higher than the length of star map.")

        # Get the star SimobjID
        brightStar = list(neighboringStarMapOnSingleSensor.SimobjID)[index]
        neighboringStar = neighboringStarMapOnSingleSensor.SimobjID[brightStar]

        # Get all star SimobjID list
        allStar = neighboringStar[:]
        allStar.append(brightStar)

        # Get the pixel positions
        allStarPosX = []
        allStarPosY = []
        for star in allStar:

            # Get the star pixel position 
            starX, starY = neighboringStarMapOnSingleSensor.RaDeclInPixel[star]

            # Transform the coordiante from DM team to camera team
            starX, starY = self.dmXY2CamXY(starX, starY)

            allStarPosX.append(starX)
            allStarPosY.append(starY)

        # Check the ccd image dimenstion
        ccdD1, ccdD2 = ccdImg.shape

        # Define the range of image
        # Get min/ max of x, y
        minX = int(min(allStarPosX))
        maxX = int(max(allStarPosX))

        minY = int(min(allStarPosY))
        maxY = int(max(allStarPosY))

        # Get the central point
        cenX = int(np.mean([minX, maxX]))
        cenY = int(np.mean([minY, maxY]))

        # Get the image dimension
        d1 = (maxY-minY) + 4*self.donutRadiusInPixel
        d2 = (maxX-minX) + 4*self.donutRadiusInPixel

        # Make d1 and d2 to be symmetric and even
        d = max(d1, d2)
        if (d%2 == 1):
            # Use d-1 instead of d+1 to avoid the boundary touch
            d = d-1

        # Compare the distances from the central point to four boundaries of ccd image
        cenYup = ccdD1 - cenY
        cenXright = ccdD2 - cenX

        # If central x or y plus d/2 will over the boundary, shift the central x, y values
        cenY = self.__shiftCenter(cenY, ccdD1, d/2)
        cenY = self.__shiftCenter(cenY, 0, d/2)

        cenX = self.__shiftCenter(cenX, ccdD2, d/2)
        cenX = self.__shiftCenter(cenX, 0, d/2)

        # Get the bright star and neighboring stas image
        offsetX = cenX-d/2
        offsetY = cenY-d/2
        singleSciNeiImg = ccdImg[int(offsetY):int(cenY+d/2), int(offsetX):int(cenX+d/2)]

        # Get the stars position in the new coordinate system
        # The final one is the bright star
        allStarPosX = np.array(allStarPosX)-offsetX
        allStarPosY = np.array(allStarPosY)-offsetY

        # Get the star magnitude
        magList = getattr(neighboringStarMapOnSingleSensor, "LSSTMag"+aFilter.upper())

        # Get the list of magnitude
        magRatio = np.array([])
        for star in allStar:
            neiMag = magList[star]
            magRatio = np.append(magRatio, neiMag)

        # Calculate the magnitude ratio
        magRatio = 1/100**((magRatio-magRatio[-1])/5.0)
        magRatio = magRatio.tolist()

        return singleSciNeiImg, allStarPosX, allStarPosY, magRatio, offsetX, offsetY

    def __shiftCenter(self, center, boundary, distance):
        """

        Shift the center if its distance to boundary is less than required.

        Arguments:
            center {[float]} -- Center point.
            boundary {[float]} -- Boundary point.
            distance {[float]} -- Required distance.

        Returns:
            [float] -- Shifted center.
        """

        # Distance between the center and boundary
        delta = boundary - center

        # Shift the center if needed
        if (abs(delta) < distance):
            center = boundary - np.sign(delta)*distance

        return center

    def simulateImg(self, imageFolderPath, defocalDis, neighboringStarMapOnSingleSensor, aFilterType, 
                    noiseRatio=0.01):
        """

        Simulate the defocal CCD images with the neighboring star map.

        Arguments:
            imageFolderPath {[str]} -- Path to image directory.
            defocalDis {[float]} -- Defocal distance in mm.
            neighboringStarMapOnSingleSensor {[dict]} -- Neighboring star map.
            aFilterType {[string]} -- Active filter type.

        Keyword Arguments:
            noiseRatio {[float]} -- The noise ratio. (default: {0.01})

        Returns:
            [float] -- Simulated intra- and extra-focal images.

        Raises:
            ValueError -- No intra-focal image files.
            ValueError -- Numbers of intra- and extra-focal image files are different.
        """

        # Generate the intra- and extra-focal ccd images
        d1, d2 = self.sensorDimList[self.sensorName]
        ccdImgIntra = np.random.random([d2, d1])*noiseRatio
        ccdImgExtra = ccdImgIntra.copy()

        # Redefine the format of defocal distance
        defocalDis = "%.2f" % defocalDis

        # Get all files in the image directory in a sorted order
        fileList = sorted(os.listdir(imageFolderPath))

        # Get the available donut files
        intraFileList = []
        extraFileList = []
        for afile in fileList:

            # Get the file name
            fileName, fileExtension = os.path.splitext(afile)

            # Split the file name for the analysis
            fileNameStr = fileName.split("_")

            # Find the file name with the correct defocal distance
            if (len(fileNameStr) == 3 and fileNameStr[1] == defocalDis):

                # Collect the file name based on the defocal type
                if (fileNameStr[-1] == "intra"):
                    intraFileList.append(afile)
                elif (fileNameStr[-1] == "extra"):
                    extraFileList.append(afile)

        # Get the number of available files
        numFile = len(intraFileList)
        if (numFile == 0):
            raise ValueError("No available donut images.")

        # Check the numbers of intra- and extra-focal images should be the same
        if (numFile != len(extraFileList)):
            raise ValueError("The numbers of intra- and extra-focal images are different.")

        # Get the magnitude of stars
        nameOfMagAttribute = "LSSTMag" + aFilterType.upper()
        starMag = getattr(neighboringStarMapOnSingleSensor, nameOfMagAttribute)

        # Based on the neighboringStarMapOnSingleSensor to reconstruct the image
        for brightStar, neighboringStar in neighboringStarMapOnSingleSensor.SimobjID.items():

            # Generate a random number
            randNum = np.random.randint(0, high=numFile)

            # Choose a random donut image from the file
            donutImageIntra = self.__getDonutImgFromFile(imageFolderPath, intraFileList[randNum])
            donutImageExtra = self.__getDonutImgFromFile(imageFolderPath, extraFileList[randNum])

            # Get the bright star magnitude
            magBS = starMag[brightStar]

            # Combine the bright star and neighboring stars. Put the bright star in the first one.
            allStars = neighboringStar[:]
            allStars.insert(0, brightStar)

            # Add the donut image
            for star in allStars:

                # Get the brigtstar pixel x, y
                starX, starY = neighboringStarMapOnSingleSensor.RaDeclInPixel[star]
                magStar = starMag[star]

                # Transform the coordiante from DM team to camera team
                starX, starY = self.dmXY2CamXY(starX, starY)

                # Ratio of magnitude between donuts (If the magnitudes of stars differs by 5,
                # the brightness differs by 100.)
                # (Magnitude difference shoulbe be >= 1.)
                magDiff = magStar-magBS
                magRatio = 1/100**(magDiff/5.0)

                # Add the donut image
                self.__addDonutImage(magRatio*donutImageIntra, starX, starY, ccdImgIntra)
                self.__addDonutImage(magRatio*donutImageExtra, starX, starY, ccdImgExtra)

        return ccdImgIntra, ccdImgExtra

    def __getDonutImgFromFile(self, imageFolderPath, fileName):
        """

        Read the donut image from the file.

        Arguments:
            imageFolderPath {[str]} -- Path to image directory.
            fileName {[str]} -- File name.

        Returns:
            [float] -- Image in numpy array.
        """

        # Get the donut image from the file by the delegation
        self.blendedImageDecorator.setImg(imageFile=os.path.join(imageFolderPath, fileName))

        return self.blendedImageDecorator.image.copy()

    def __addDonutImage(self, donutImage, starX, starY, ccdImg):
        """

        Add the donut image to simulated CCD image frame.

        Arguments:
            donutImage {[float]} -- Image in numpy array.
            starX {[float]} -- Star position in pixel x.
            starY {[float]} -- Star position in pixel y.
            ccdImg {[float]} -- CCD image in numpy array.
        """

        # Get the dimension of donut image
        d1, d2 = donutImage.shape

        # Get the interger of position to use as the index
        y = int(starY)
        x = int(starX)

        # Add the donut image on the CCD image
        ccdImg[y-int(d1/2):y-int(d1/2)+d1, x-int(d2/2):x-int(d2/2)+d2] += donutImage

def expandDetectorName(abbrevName):
    """Convert a detector name of the form Rxy_Sxy[_Ci] to canonical form: R:x,y S:x,y[,c]
    C0 -> A, C1 -> B

    This is copied from lsst.obs.lsstSim:
    https://github.com/lsst/obs_lsstSim/blob/master/bin.src/makeLsstCameraRepository.py
    """
    m = re.match(r"R(\d)(\d)_S(\d)(\d)(?:_C([0,1]))?$", abbrevName)
    if m is None:
        raise RuntimeError("Cannot parse abbreviated name %r" % (abbrevName,))
    fullName = "R:%s,%s S:%s,%s" % tuple(m.groups()[0:4])
    subSensor = m.groups()[4]
    if subSensor is not None:
        fullName = fullName + "," + {"0": "A", "1": "B"}[subSensor]
    return fullName

def abbrevDectectorName(canonicalForm):
    """
    
    Convert a canonical name to abbreviate name (R:x,y S:x,y[,c] --> Rxy_Sxy[_Ci]).
    
    Arguments:
        canonicalForm {[str]} -- Detector canonical name.
    
    Returns:
        [str] -- Abbreviated name.
    
    Raises:
        RuntimeError -- Input does not match the canonical form of detector name.
    """

    # Use the regular expression to analyze the input name
    m = re.match(r"R:(\d),(\d) S:(\d),(\d)(?:,([A,B]))?$", canonicalForm)

    # Raise error if the input does not match the form of regular expression
    if m is None:
        raise RuntimeError("Cannot parse canonical name %r" % (canonicalForm,))

    # Generate the abbreviated name
    abbrevName = "R%s%s_S%s%s" % tuple(m.groups()[0:4])

    # Check the sensor is wavefront sensor or not
    subSensor = m.groups()[4]
    if (subSensor is not None):
        # Label the WFS sensor
        abbrevName = abbrevName + "_C" + {"A": "0", "B": "1"}[subSensor]

    # Return the abbreviate name
    return abbrevName


if __name__ == "__main__":
    pass
