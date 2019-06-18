import re
import numpy as np

from lsst.ts.wep.ButlerWrapper import ButlerWrapper
from lsst.ts.wep.DefocalImage import DefocalImage
from lsst.ts.wep.DonutImage import DonutImage
from lsst.ts.wep.Utility import abbrevDectectorName, searchDonutPos, \
    DefocalType, ImageType


class WepController(object):

    CORNER_WFS_LIST = ["R:0,0 S:2,2,A", "R:0,0 S:2,2,B", "R:0,4 S:2,0,A",
                       "R:0,4 S:2,0,B", "R:4,0 S:0,2,A", "R:4,0 S:0,2,B",
                       "R:4,4 S:0,0,A", "R:4,4 S:0,0,B"]

    def __init__(self, dataCollector, isrWrapper, sourSelc, sourProc, wfEsti):
        """Initialize the wavefront estimation pipeline (WEP) controller class.

        Parameters
        ----------
        dataCollector : CamDataCollector
            Camera data collector.
        isrWrapper : CamIsrWrapper
            Instrument signature removal (ISR) wrapper.
        sourSelc : SourceSelector
            Source selector.
        sourProc : SourceProcessor
            Source processor.
        wfEsti : WfEstimator
            Wavefront estimator.
        """

        # Data collector to use DM ingest command line task
        self.dataCollector = dataCollector

        # ISR wrapper to use DM ISR command line task
        self.isrWrapper = isrWrapper

        # Source selector to query the bright star catalog
        self.sourSelc = sourSelc

        # Source processor to get the donut image or do the deblending
        self.sourProc = sourProc

        # Wavefront estimator to get zk
        self.wfEsti = wfEsti

        # Butler wrapper to use DM data butler
        self.butlerWrapper = None

    def getDataCollector(self):
        """Get the attribute of data collector.

        Returns
        -------
        CamDataCollector
            Data collector.
        """

        return self.dataCollector

    def getIsrWrapper(self):
        """Get the attribute of ISR wraper.

        ISR: Instrument signature removal.

        Returns
        -------
        CamIsrWrapper
            ISR wraper.
        """

        return self.isrWrapper

    def getSourSelc(self):
        """Get the attribute of source selector.

        Returns
        -------
        SourceSelector
            Source selector.
        """

        return self.sourSelc

    def getSourProc(self):
        """Get the attribute of source processor.

        Returns
        -------
        SourceProcessor
            Source processor.
        """

        return self.sourProc

    def getWfEsti(self):
        """Get the attribute of wavefront estimator.

        Returns
        -------
        WfEstimator
            Wavefront estimator.
        """

        return self.wfEsti

    def getButlerWrapper(self):
        """Get the attribute of butler wrapper.

        Returns
        -------
        ButlerWrapper
            Butler wrapper.
        """

        return self.butlerWrapper

    def setPostIsrCcdInputs(self, inputs):
        """Set inputs of post instrument signature removal (ISR) CCD images.

        Parameters
        ----------
        inputs : RepositoryArgs, dict, or str
            Can be a single item or a list. Provides arguments to load an
            existing repository (or repositories). String is assumed to be a
            URI and is used as the cfgRoot (URI to the location of the cfg
            file). (Local file system URI does not have to start with
            'file://' and in this way can be a relative path). The
            'RepositoryArgs' class can be used to provide more parameters with
            which to initialize a repository (such as 'mapper', 'mapperArgs',
            'tags', etc. See the 'RepositoryArgs' documentation for more
            details). A dict may be used as shorthand for a 'RepositoryArgs'
            class instance. The dict keys must match parameters to the
            'RepositoryArgs.__init__' function.
        """

        self.butlerWrapper = ButlerWrapper(inputs)

    def getPostIsrImgMapByPistonDefocal(self, sensorNameList, obsIdList):
        """Get the post ISR image map that the defocal images are by the
        pistion motion.

        Parameters
        ----------
        sensorNameList : list
            List of sensor name.
        obsIdList : list
            Observation Id list in [intraObsId, extraObsId].

        Returns
        -------
        dict
            Post-ISR image map. The dictionary key is the sensor name. The
            dictionary item is the defocal image on the camera coordinate.
            (type: DefocalImage).
        """

        return self._getImgMapByPistonDefocal(sensorNameList, obsIdList,
                                              ImageType.Amp)

    def _getImgMapByPistonDefocal(self, sensorNameList, obsIdList, imageType):
        """Get the image map that the defocal images are by the pistion motion.

        Parameters
        ----------
        sensorNameList : list
            List of sensor name.
        obsIdList : list
            Observation Id list in [intraObsId, extraObsId].
        imageType : ImageType
            Image type.

        Returns
        -------
        dict
            Image map. The dictionary key is the sensor name. The dictionary
            item is the defocal image on the camera coordinate. (type:
            DefocalImage).

        Raises
        ------
        ValueError
            The image type is not supported.
        """

        # Get the waveront image map
        wfsImgMap = dict()
        for sensorName in sensorNameList:

            # Get the sensor name information
            raft, sensor = self._getSensorInfo(sensorName)[0:2]

            # The intra/ extra defocal images are decided by obsId
            imgList = []
            for visit in obsIdList:

                # Get the exposure image in ndarray
                if (imageType == ImageType.Amp):
                    exp = self.butlerWrapper.getPostIsrCcd(
                        int(visit), raft, sensor)
                elif (imageType == ImageType.Eimg):
                    exp = self.butlerWrapper.getEimage(
                        int(visit), raft, sensor)
                else:
                    raise ValueError("The %s is not supported." % imageType)

                img = self.butlerWrapper.getImageData(exp)

                # Transform the image in DM coordinate to camera coordinate.
                camImg = self._transImgDmCoorToCamCoor(img)

                # Collect the image
                imgList.append(camImg)

            wfsImgMap[sensorName] = DefocalImage(intraImg=imgList[0],
                                                 extraImg=imgList[1])

        return wfsImgMap

    def _getSensorInfo(self, sensorName):
        """Get the sensor information.

        Parameters
        ----------
        sensorName : str
            Sensor name (e.g. "R:2,2 S:1,1" or "R:0,0 S:2,2,A")

        Returns
        -------
        str
            Raft.
        str
            Sensor.
        str
            Channel.
        """

        raft = sensor = channel = None

        # Use the regular expression to analyze the input name
        m = re.match(r"R:(\d,\d) S:(\d,\d)(?:,([A,B]))?$", sensorName)
        if (m is not None):
            raft, sensor, channel = m.groups()[0:3]

        # This is for the phosim mapper use.
        # For example, raft is "R22" and sensor is "S11".
        raftAbbName = "R" + raft[0] + raft[-1]
        sensorAbbName = "S" + sensor[0] + sensor[-1]

        return raftAbbName, sensorAbbName, channel

    def _transImgDmCoorToCamCoor(self, dmImg):
        """Transfrom the image in DM coordinate to camera coordinate.

        The geometry transformation is defined in LSE-349.

        Parameters
        ----------
        dmImg : numpy.ndarray
            Image in DM coordinate.

        Returns
        -------
        numpy.ndarray
            Image is camera coordinate.
        """

        # For the PhoSim mapper, it is the transpose (x, y) to (y, x).
        camImg = dmImg.T

        return camImg

    def getEimgMapByPistonDefocal(self, sensorNameList, obsIdList):
        """Get the eimage map that the defocal images are by the pistion
        motion.

        Parameters
        ----------
        sensorNameList : list
            List of sensor name.
        obsIdList : list
            Observation Id list in [intraObsId, extraObsId].

        Returns
        -------
        dict
            Eimage map. The dictionary key is the sensor name. The dictionary
            item is the defocal image on the camera coordinate. (type:
            DefocalImage).
        """

        return self._getImgMapByPistonDefocal(sensorNameList, obsIdList,
                                              ImageType.Eimg)

    def getPostIsrImgMapOnCornerWfs(self, sensorNameList, obsId):
        """Get the post ISR image map of corner wavefront sensors.

        Parameters
        ----------
        sensorNameList : list
            List of sensor name.
        obsId : int
            Observation Id.

        Returns
        -------
        dict
            Post-ISR image map. The dictionary key is the sensor name. The
            dictionary item is the defocal image on the camera coordinate.
            (type: DefocalImage).
        """
        pass

    def getDonutMap(self, neighborStarMap, wfsImgMap, filterType,
                    doDeblending=False):
        """Get the donut map on each wavefront sensor (WFS).

        Parameters
        ----------
        neighborStarMap : dict
            Information of neighboring stars and candidate stars with the name
            of sensor as a dictionary.
        wfsImgMap : dict
            Post-ISR image map. The dictionary key is the sensor name. The
            dictionary item is the defocal image on the camera coordinate.
            (type: DefocalImage).
        filterType : FilterType
            Filter type.
        doDeblending : bool, optional
            Do the deblending or not. If False, only consider the single donut
            based on the bright star catalog.(the default is False.)

        Returns
        -------
        dict
            Donut image map. The dictionary key is the sensor name. The
            dictionary item is the list of donut image (type:
            list[DonutImage]).
        """

        donutMap = dict()
        for sensorName, nbrStar in neighborStarMap.items():

            # Get the abbraviated sensor name
            abbrevName = abbrevDectectorName(sensorName)

            # Configure the source processor
            self.sourProc.config(sensorName=abbrevName)

            # Get the defocal images: [intra, extra]
            defocalImgList = [wfsImgMap[sensorName].getIntraImg(),
                              wfsImgMap[sensorName].getExtraImg()]

            # Get the bright star id list on specific sensor
            brightStarIdList = list(nbrStar.getId())
            for starIdIdx in range(len(brightStarIdList)):

                # Get the single star map
                for jj in range(len(defocalImgList)):

                    ccdImg = defocalImgList[jj]

                    # Get the segment of image
                    if (ccdImg is not None):
                        singleSciNeiImg, allStarPosX, allStarPosY, magRatio, \
                            offsetX, offsetY = \
                            self.sourProc.getSingleTargetImage(
                                ccdImg, nbrStar, starIdIdx, filterType)

                        # Only consider the single donut if no deblending
                        if (not doDeblending) and (len(magRatio) != 1):
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
                            imgDeblend, realcx, realcy = \
                                self.sourProc.doDeblending(
                                    singleSciNeiImg, allStarPosX, allStarPosY,
                                    magRatio)
                            # Update the magnitude ratio
                            magRatio = [1]

                        else:
                            continue

                        # Extract the image
                        if (len(magRatio) == 1):
                            sizeInPix = self.wfEsti.getSizeInPix()
                            x0 = np.floor(realcx - sizeInPix / 2).astype("int")
                            y0 = np.floor(realcy - sizeInPix / 2).astype("int")
                            imgDeblend = imgDeblend[y0:y0 + sizeInPix,
                                                    x0:x0 + sizeInPix]

                        # Rotate the image if the sensor is the corner
                        # wavefront sensor
                        if sensorName in self.CORNER_WFS_LIST:

                            # Get the Euler angle
                            eulerZangle = round(self.sourProc.getEulerZinDeg(
                                abbrevName))

                            # Change the sign if the angle < 0
                            while (eulerZangle < 0):
                                eulerZangle += 360

                            # Do the rotation of matrix
                            numOfRot90 = eulerZangle // 90
                            imgDeblend = np.flipud(
                                np.rot90(np.flipud(imgDeblend), numOfRot90))

                        # Put the deblended image into the donut map
                        if sensorName not in donutMap.keys():
                            donutMap[sensorName] = []

                        # Check the donut exists in the list or not
                        starId = brightStarIdList[starIdIdx]
                        donutIndex = self._searchDonutListId(
                            donutMap[sensorName], starId)

                        # Create the donut object and put into the list if it
                        # is needed
                        if (donutIndex < 0):

                            # Calculate the field X, Y
                            pixelX = realcx + offsetX
                            pixelY = realcy + offsetY
                            fieldX, fieldY = self.sourProc.camXYtoFieldXY(
                                pixelX, pixelY)

                            # Instantiate the DonutImage class
                            donutImg = DonutImage(starId, pixelX, pixelY,
                                                  fieldX, fieldY)
                            donutMap[sensorName].append(donutImg)

                            # Search for the donut index again
                            donutIndex = self._searchDonutListId(
                                donutMap[sensorName], starId)

                        # Get the donut image list
                        donutList = donutMap[sensorName]

                        # Take the absolute value for images, which might
                        # contain the negative value after the ISR correction.
                        # This happens for the amplifier images.
                        imgDeblend = np.abs(imgDeblend)

                        # Set the intra focal image
                        if (jj == 0):
                            donutList[donutIndex].setImg(intraImg=imgDeblend)
                        # Set the extra focal image
                        elif (jj == 1):
                            donutList[donutIndex].setImg(extraImg=imgDeblend)

        return donutMap

    def _searchDonutListId(self, donutList, starId):
        """Search the bright star ID in the donut list.

        Parameters
        ----------
        donutList : list
            List of DonutImage object.
        starId : int
            Star Id.

        Returns
        -------
        int
             Index of donut image object with the specific starId.
        """

        index = -1
        for ii in range(len(donutList)):
            if (donutList[ii].getStarId() == int(starId)):
                index = ii
                break

        return index

    def calcWfErr(self, donutMap):
        """Calculate the wavefront error in annular Zernike polynomials
        (z4-z22).

        Parameters
        ----------
        donutMap : dict
            Donut image map. The dictionary key is the sensor name. The
            dictionary item is the donut image (type: DonutImage).

        Returns
        -------
        dict
            Donut image map with calculated wavefront error.
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
                intraFieldXY = intraDonut.getFieldPos()
                extraFieldXY = extraDonut.getFieldPos()

                # Get the defocal images
                intraImg = intraDonut.getIntraImg()
                extraImg = extraDonut.getExtraImg()

                # Calculate the wavefront error
                zer4UpNm = self._calcSglWfErr(intraImg, extraImg, intraFieldXY,
                                              extraFieldXY)

                # Put the value to the donut image
                intraDonut.setWfErr(zer4UpNm)
                extraDonut.setWfErr(zer4UpNm)

        # Intentionally to expose this return value to show the input,
        # donutMap, has been modified.
        return donutMap

    def _calcSglWfErr(self, intraImg, extraImg, intraFieldXY, extraFieldXY):
        """Calculate the wavefront error in annular Zernike polynomials
        (z4-z22) for single donut.

        Parameters
        ----------
        intraImg : numpy.ndarray
            Intra-focal donut image.
        extraImg : numpy.ndarray
            Extra-focal donut image.
        intraFieldXY : tuple
            Field x, y in degree of intra-focal donut image.
        extraFieldXY : tuple
            Field x, y in degree of extra-focal donut image.

        Returns
        -------
        numpy.ndarray
            Coefficients of Zernike polynomials (z4 - z22) in nm.
        """

        # Set the images
        self.wfEsti.setImg(intraFieldXY, DefocalType.Intra, image=intraImg)
        self.wfEsti.setImg(extraFieldXY, DefocalType.Extra, image=extraImg)

        # Reset the wavefront estimator
        self.wfEsti.reset()

        # Calculate the wavefront error
        zer4UpNm = self.wfEsti.calWfsErr()

        return zer4UpNm

    def calcAvgWfErrOnSglCcd(self, donutList):
        """Calculate the average of wavefront error on single CCD.

        CCD: Charge-coupled device.

        Parameters
        ----------
        donutList : list
            List of donut object (type: DonutImage).

        Returns
        -------
        numpy.ndarray
            Average of wavefront error in nm.
        """

        # Calculate the weighting of donut image
        wgtRatio = self._calcWeiRatio(donutList)

        # Calculate the mean wavefront error
        avgErr = 0
        for ii in range(len(donutList)):

            donut = donutList[ii]
            zer4UpNmArr = donut.getWfErr()

            # Assign the zer4UpNm
            if (len(zer4UpNmArr) == 0):
                zer4UpNm = 0
            else:
                zer4UpNm = zer4UpNmArr

            avgErr = avgErr + wgtRatio[ii] * zer4UpNm

        return avgErr

    def _calcWeiRatio(self, donutList):
        """Calculate the weighting ratio of donut in the list.

        Parameters
        ----------
        donutList : list
            List of donut object (type: DonutImage).

        Returns
        -------
        numpy.ndarray
            Array of weighting ratio of donuts to do the average of wavefront
            error.
        """

        # Weighting of donut image. Use the simple average at this moment.
        # Need to consider the S/N and other factors in the future.

        # Check the available zk and give the ratio
        wgtRatio = []
        for donut in donutList:
            if (len(donut.getWfErr()) == 0):
                wgtRatio.append(0)
            else:
                wgtRatio.append(1)

        # Do the normalization
        wgtRatioArr = np.array(wgtRatio)
        normalizedwgtRatioArr = wgtRatioArr / np.sum(wgtRatioArr)

        return normalizedwgtRatioArr

    def genMasterDonut(self, donutMap, zcCol=np.zeros(22)):
        """Generate the master donut map.

        Parameters
        ----------
        donutMap : dict
            Donut image map. The dictionary key is the sensor name. The
            dictionary item is the donut image (type: DonutImage).
        zcCol : numpy.ndarray, optional
            Coefficients of wavefront (z1-z22) in nm. (the default is
            np.zeros(22).)

        Returns
        -------
        dict
            Master donut image map. The dictionary key is the sensor name. The
            dictionary item is the master donut image (type: DonutImage).
        """

        masterDonutMap = dict()
        for sensorName, donutList in donutMap.items():

            # Get the master donut on single CCD
            masterDonut = self._genMasterImgOnSglCcd(donutList, zcCol=zcCol)

            # Put the master donut to donut map
            masterDonutMap[sensorName] = [masterDonut]

        return masterDonutMap

    def _genMasterImgOnSglCcd(self, donutList, zcCol):
        """Generate the master donut image on single CCD.

        CCD: Charge-coupled device.

        Parameters
        ----------
        sensorName : str
            Canonical sensor name (e.g. "R:2,2 S:1,1").
        donutList : list
            List of donut object (type: DonutImage).
        zcCol : numpy.ndarray
            Coefficients of wavefront (z1-z22) in nm.

        Returns
        -------
        DonutImage
            Master donut.
        """

        intraProjImgList = []
        extraProjImgList = []
        for donut in donutList:

            # Get the field x, y
            fieldXY = donut.getFieldPos()

            # Set the image
            intraImg = donut.getIntraImg()
            if (intraImg is not None):

                # Get the projected image
                projImg = self._getProjImg(fieldXY, DefocalType.Intra,
                                           intraImg, zcCol)

                # Collect the projected donut
                intraProjImgList.append(projImg)

            extraImg = donut.getExtraImg()
            if (extraImg is not None):

                # Get the projected image
                projImg = self._getProjImg(fieldXY, DefocalType.Extra,
                                           extraImg, zcCol)

                # Collect the projected donut
                extraProjImgList.append(projImg)

        # Generate the master donut
        stackIntraImg = self._stackImg(intraProjImgList)
        stackExtraImg = self._stackImg(extraProjImgList)

        # Put the master donut to donut map
        pixelX, pixelY = searchDonutPos(stackIntraImg)
        masterDonut = DonutImage(0, pixelX, pixelY, 0, 0,
                                 intraImg=stackIntraImg, extraImg=stackExtraImg)

        return masterDonut

    def _getProjImg(self, fieldXY, defocalType, defocalImg, zcCol):
        """Get the projected image on the pupil.

        Parameters
        ----------
        fieldXY : tuple
            Position of donut on the focal plane in the degree for intra- and
            extra-focal images (field X, field Y).
        defocalType : enum 'DefocalType'
            Defocal type of image.
        defocalImg : numpy.ndarray
            Donut defocal image.
        zcCol : numpy.ndarray
            Coefficients of wavefront (z1-z22).

        Returns
        -------
        numpy.ndarray
            Projected image.
        """

        # Set the image
        self.wfEsti.setImg(fieldXY, defocalType, image=defocalImg)

        # Get the image in the type of CompensationImageDecorator
        if (defocalType == DefocalType.Intra):
            img = self.wfEsti.getIntraImg()
        elif (defocalType == DefocalType.Extra):
            img = self.wfEsti.getExtraImg()

        # Get the distortion correction (offaxis)
        algo = self.wfEsti.getAlgo()
        offAxisCorrOrder = algo.getOffAxisPolyOrder()

        inst = self.wfEsti.getInst()
        img.setOffAxisCorr(inst, offAxisCorrOrder)

        # Do the image cocenter
        img.imageCoCenter(inst)

        # Do the compensation/ projection
        img.compensate(inst, algo, zcCol, self.wfEsti.getOptModel())

        # Return the projected image
        return img.getImg()

    def _stackImg(self, imgList):
        """Stack the images.

        Parameters
        ----------
        imgList : list
            List of image.

        Returns
        -------
        numpy.ndarray
            Stacked image.
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
                cy = int(dy / 2)
                cx = int(dx / 2)
                stackImg += img[cy - deltaY:cy + deltaY,
                                cx - deltaX:cx + deltaX]

        return stackImg


if __name__ == "__main__":
    pass
