import os
import numpy as np

from lsst.ts.wep.deblend.BlendedImageDecorator import BlendedImageDecorator
from lsst.ts.wep.Utility import readPhoSimSettingData


class SourceProcessor(object):

    # Donut radius is 63 pixel if the defocal distance is 1.5 mm
    STAR_RADIUS_IN_PIXEL = 63

    # 1 pixel = 0.2 arcsec
    PIXEL_TO_ARCSEC = 0.2

    # 1 pixel = 10 um
    PIXEL_TO_UM = 10

    PHOSIM_FOCALPLANE = "focalplanelayout.txt"

    def __init__(self, folderPath2FocalPlane=None):
        """Initialize the SourceProcessor class.

        Parameters
        ----------
        folderPath2FocalPlane : str, optional
            Directory of focal plane file. (the default is None.)
        """

        self.sensorName = ""
        self.blendedImageDecorator = BlendedImageDecorator()

        self.sensorFocaPlaneInDeg = dict()
        self.sensorFocaPlaneInUm = dict()
        self.sensorDimList = dict()
        self.sensorEulerRot = dict()

        if (folderPath2FocalPlane is not None):
            self._readFocalPlane(folderPath2FocalPlane)

    def _readFocalPlane(self, folderPath):
        """Read the focal plane data used in PhoSim to get the ccd dimension
        and fieldXY in chip center.

        Parameters
        ----------
        folderPath : str
            Directory of focal plane file.
        """

        # Read the focal plane data by the delegation
        ccdData = readPhoSimSettingData(folderPath, self.PHOSIM_FOCALPLANE,
                                        "fieldCenter")

        # Collect the focal plane data
        sensorFocaPlaneInDeg = dict()
        sensorFocaPlaneInUm = dict()
        sensorDimList = dict()
        for sensorName, data in ccdData.items():

            # Consider the x-translation in corner wavefront sensors
            self._shiftCenterWfs(sensorName, data)

            # Change the unit from um to degree
            xInUm = float(data[0])
            yInUm = float(data[1])
            pixelSizeInUm = float(data[2])
            sizeXinPixel = int(data[3])
            sizeYinPixel = int(data[4])

            # 1 degree = 3600 arcsec
            fieldX = xInUm / pixelSizeInUm * self.PIXEL_TO_ARCSEC / 3600
            fieldY = yInUm / pixelSizeInUm * self.PIXEL_TO_ARCSEC / 3600

            # Get the data
            sensorFocaPlaneInDeg[sensorName] = (fieldX, fieldY)
            sensorFocaPlaneInUm[sensorName] = (xInUm, yInUm)
            sensorDimList[sensorName] = (sizeXinPixel, sizeYinPixel)

        # Assign the values
        self.sensorDimList = sensorDimList
        self.sensorFocaPlaneInDeg = sensorFocaPlaneInDeg
        self.sensorFocaPlaneInUm = sensorFocaPlaneInUm
        self.sensorEulerRot = readPhoSimSettingData(
                                folderPath, self.PHOSIM_FOCALPLANE, "eulerRot")

    def _shiftCenterWfs(self, sensorName, focalPlaneData):
        """Shift the fieldXY of center of wavefront sensors.

        The input data is the center of combined chips (C0+C1). The input data
        will be updated directly.

        Parameters
        ----------
        sensorName : str
            Abbreviated sensor name.
        focalPlaneData : list
            Data of focal plane: [x position (microns), y position (microns),
            pixel size (microns), number of x pixels, number of y pixels].
        """

        # The layout is shown in the following:

        # R04_S20              R44_S00
        # --------           -----------       /\ +y
        # |  C0  |           |    |    |        |
        # |------|           | C1 | C0 |        |
        # |  C1  |           |    |    |        |
        # --------           -----------        -----> +x

        # R00_S22              R40_S02
        # -----------          --------
        # |    |    |          |  C1  |
        # | C0 | C1 |          |------|
        # |    |    |          |  C0  |
        # -----------          --------

        xInUm = float(focalPlaneData[0])
        yInUm = float(focalPlaneData[1])
        pixelSizeInUm = float(focalPlaneData[2])
        sizeXinPixel = float(focalPlaneData[3])

        # Consider the x-translation in corner wavefront sensors
        tempX = None
        tempY = None

        if sensorName in ("R44_S00_C0", "R00_S22_C1"):
            # Shift center to +x direction
            tempX = xInUm + sizeXinPixel / 2 * pixelSizeInUm
        elif sensorName in ("R44_S00_C1", "R00_S22_C0"):
            # Shift center to -x direction
            tempX = xInUm - sizeXinPixel / 2 * pixelSizeInUm
        elif sensorName in ("R04_S20_C1", "R40_S02_C0"):
            # Shift center to -y direction
            tempY = yInUm - sizeXinPixel / 2 * pixelSizeInUm
        elif sensorName in ("R04_S20_C0", "R40_S02_C1"):
            # Shift center to +y direction
            tempY = yInUm + sizeXinPixel / 2 * pixelSizeInUm

        # Replace the value by the shifted one
        if (tempX is not None):
            focalPlaneData[0] = str(tempX)
        elif (tempY is not None):
            focalPlaneData[1] = str(tempY)

    def config(self, sensorName=None, folderPath2FocalPlane=None):
        """Do the configuration.

        Parameters
        ----------
        sensorName : str, optional
            Abbreviated sensor name. (the default is None.)
        folderPath2FocalPlane : str, optional
            Directory of focal plane file. (the default is None.)
        """

        # Give the sensor name
        if (sensorName is not None):
            self.sensorName = sensorName

        # Read the focal plane data
        if (folderPath2FocalPlane is not None):
            self._readFocalPlane(folderPath2FocalPlane)

    def getEulerZinDeg(self, sensorName):
        """Get the Euler Z angle of sensor in degree.

        Parameters
        ----------
        sensorName : str
            Abbreviated sensor name.

        Returns
        -------
        float
            Euler Z angle in degree.
        """

        return float(self.sensorEulerRot[sensorName][0])

    def camXYtoFieldXY(self, pixelX, pixelY):
        """Get the field X, Y from the pixel x, y position on CCD.

        Parameters
        ----------
        pixelX : float
            Pixel x on camera coordinate.
        pixelY : float
            Pixel y on camera coordinate.

        Returns
        -------
        float
            Field x in degree.
        float
            Field y in degree.
        """

        # The wavefront sensors will do the counter-clockwise rotation as the
        # following based on the euler angle:

        # R04_S20              R44_S00
        # O-------           -----O----O       /\ +y
        # |  C0  |           |    |    |        |
        # O------|           | C1 | C0 |        |
        # |  C1  |           |    |    |        |
        # --------           -----------        O----> +x

        # R00_S22              R40_S02
        # -----------          --------
        # |    |    |          |  C1  |
        # | C0 | C1 |          |------O
        # |    |    |          |  C0  |
        # O----O-----          -------O

        # Get the field X, Y of sensor's center
        fieldXc, fieldYc = self.sensorFocaPlaneInDeg[self.sensorName]

        # Get the center pixel position
        pixelXc, pixelYc = self.sensorDimList[self.sensorName]
        pixelXc = pixelXc / 2
        pixelYc = pixelYc / 2

        # Calculate the delta x and y in degree
        # 1 degree = 3600 arcsec
        deltaX = (pixelX - pixelXc) * self.PIXEL_TO_ARCSEC / 3600.0
        deltaY = (pixelY - pixelYc) * self.PIXEL_TO_ARCSEC / 3600.0

        # Calculate the transformed coordinate in degree.
        fieldX, fieldY = self._rotCam2FocalPlane(
                            self.sensorName, fieldXc, fieldYc, deltaX, deltaY)

        return fieldX, fieldY

    def _rotCam2FocalPlane(self, sensorName, centerX, centerY, deltaX, deltaY, 
                           clockWise=False):
        """Do the rotation from camera coordinate to focal plane coordinate or
        vice versa.

        Parameters
        ----------
        sensorName : str
            Abbreviated sensor name.
        centerX : float
            CCD center x.
        centerY : float
            CCD center y.
        deltaX : float
            Delta x from the CCD's center.
        deltaY : float
            Delta y from the CCD's center.
        clockWise : bool, optional
            Rotation direction (True: clockwise, False: counter-clockwise).
            (the default is False.)

        Returns
        -------
        float
            Transformed x position.
        float
            Transformed y position.
        """

        # Get the euler angle in z direction (only consider the z rotatioin at
        # this moment)
        eulerZ = round(self.getEulerZinDeg(sensorName))
        eulerZinRad = np.deg2rad(eulerZ)

        # Counter-clockwise or clockwise rotation
        if (clockWise):
            eulerZinRad = -eulerZinRad

        # Calculate the new x, y by the rotation. This is important for
        # wavefront sensor.
        newX = centerX + np.cos(eulerZinRad) * deltaX - \
               np.sin(eulerZinRad) * deltaY
        newY = centerY + np.sin(eulerZinRad) * deltaX + \
               np.cos(eulerZinRad)*deltaY

        return newX, newY

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
        deltaX = (xInUm-xc) / self.PIXEL_TO_UM
        deltaY = (yInUm-yc) / self.PIXEL_TO_UM

        # Calculate the transformed coordinate
        pixelX, pixelY = self._rotCam2FocalPlane(self.sensorName, pixelXc, pixelYc, deltaX, 
                                                    deltaY, clockWise=True)

        return pixelX, pixelY

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
        d1 = (maxY-minY) + 4 * self.STAR_RADIUS_IN_PIXEL
        d2 = (maxX-minX) + 4 * self.STAR_RADIUS_IN_PIXEL

        # Make d1 and d2 to be symmetric and even
        d = max(d1, d2)
        if (d%2 == 1):
            # Use d-1 instead of d+1 to avoid the boundary touch
            d = d-1

        # Compare the distances from the central point to four boundaries of ccd image
        cenYup = ccdD1 - cenY
        cenXright = ccdD2 - cenX

        # If central x or y plus d/2 will over the boundary, shift the central x, y values
        cenY = self._shiftCenter(cenY, ccdD1, d/2)
        cenY = self._shiftCenter(cenY, 0, d/2)

        cenX = self._shiftCenter(cenX, ccdD2, d/2)
        cenX = self._shiftCenter(cenX, 0, d/2)

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

    def _shiftCenter(self, center, boundary, distance):
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
            donutImageIntra = self._getDonutImgFromFile(imageFolderPath, intraFileList[randNum])
            donutImageExtra = self._getDonutImgFromFile(imageFolderPath, extraFileList[randNum])

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
                self._addDonutImage(magRatio*donutImageIntra, starX, starY, ccdImgIntra)
                self._addDonutImage(magRatio*donutImageExtra, starX, starY, ccdImgExtra)

        return ccdImgIntra, ccdImgExtra

    def _getDonutImgFromFile(self, imageFolderPath, fileName):
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

    def _addDonutImage(self, donutImage, starX, starY, ccdImg):
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


if __name__ == "__main__":
    pass
