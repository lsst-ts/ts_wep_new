import os

from lsst.ts.wep.Utility import getModulePath, getConfigDir, FilterType, \
    BscDbType
from lsst.ts.wep.CamDataCollector import CamDataCollector
from lsst.ts.wep.CamIsrWrapper import CamIsrWrapper
from lsst.ts.wep.SourceProcessor import SourceProcessor
from lsst.ts.wep.SourceSelector import SourceSelector
from lsst.ts.wep.WfEstimator import WfEstimator
from lsst.ts.wep.WepController import WepController
from lsst.ts.wep.ctrlIntf.SensorWavefrontData import SensorWavefrontData


class WEPCalculation(object):
    """Base class for converting the wavefront images into wavefront errors.

    There will be different implementations of this for different
    types of CCDs (normal, full array mode, comcam, cmos, shwfs).
    """

    def __init__(self, astWcsSol, camType, isrDir):
        """Construct an WEP calculation object.

        Parameters
        ----------
        astWcsSol : AstWcsSol
            AST world coordinate system (WCS) solution.
        camType : CamType
            Camera type.
        isrDir : str
            Instrument signature remocal (ISR) directory. This directory will
            have the input and output that the data butler needs.
        """

        super().__init__()

        # This attribute is just a stakeholder here since there is no detail of
        # AST WCS solution yet
        self.astWcsSol = astWcsSol

        # ISR directory that the data butler uses
        self.isrDir = isrDir

        # Number of processors for WEP to use
        self.numOfProc = 1

        # Boresight infomation
        self.raInDeg = 0.0
        self.decInDeg = 0.0

        # Camera rotation angle
        self.rotAngInDeg = 0.0

        # Sky information file fro the temporary use
        self.skyFile = ""

        # Remove these two in a latter time
        # self.currentFilter = FilterType.REF
        # self.calibsDir = ""

        # WEP controller
        self.wepCntlr = self._configWepController(camType, isrDir)

    def _configWepController(self, camType, isrDir):
        """Configure the WEP controller.

        WEP: wavefront estimation pipeline.

        Parameters
        ----------
        camType : CamType
            Camera type.
        isrDir : str
            Instrument signature remocal (ISR) directory. This directory will
            have the input and output that the data butler needs.

        Returns
        -------
        WepController
            Configured WEP controller.
        """

        dataCollector = CamDataCollector(isrDir)
        isrWrapper = CamIsrWrapper(isrDir)
        sourSelc = self._configSourceSelector(camType)
        sourProc = self._configSourceProcessor()
        wfsEsti = self._configWfEstimator()
        wepCntlr = WepController(dataCollector, isrWrapper, sourSelc,
                                 sourProc, wfsEsti)

        return wepCntlr

    def _configSourceSelector(self, camType):
        """Configue the source selector.

        Returns
        -------
        SourceSelector
            Configured source selector.
        """

        sourSelc = SourceSelector(camType, BscDbType.LocalDbForStarFile)
        dbAdress = os.path.join(getModulePath(), "tests", "testData",
                                "bsc.db3")
        sourSelc.connect(dbAdress)

        return sourSelc

    def _configSourceProcessor(self):
        """Configure the source processor.
        Returns
        -------
        SourceProcessor
            Configured source processor.
        """

        folderPath2FocalPlane = os.path.join(getModulePath(), "tests",
                                             "testData")
        sourProc = SourceProcessor(folderPath2FocalPlane=folderPath2FocalPlane)

        return sourProc

    def _configWfEstimator(self):
        """Configure the wavefront estimator.
        Returns
        -------
        WfEstimator
            Configured wavefront estimator.
        """

        instruFolderPath = os.path.join(self._getConfigDataPath(), "cwfs",
                                        "instruData")
        algoFolderPath = os.path.join(self._getConfigDataPath(), "cwfs",
                                      "algo")
        wfsEsti = WfEstimator(instruFolderPath, algoFolderPath)

        # Use the comcam to calculate the LSST central raft image
        wfsEsti.config(solver="exp", instName="comcam",
                       opticalModel="offAxis",
                       defocalDisInMm=self.DEFOCAL_DIS_IN_MM,
                       sizeInPix=self.DONUT_IMG_SIZE_IN_PIXEL, debugLevel=0)

        return wfsEsti

    def disconnect(self):
        """Disconnect the database."""

        self.wepCntlr.getSourSelc().disconnect()

    def getIsrDir(self):
        """Get the instrument signature removal (ISR) directory.

        This directory will have the input and output that the data butler
        needs.

        Returns
        -------
        str
            ISR directory.
        """

        return self.isrDir

    def setSkyFile(self, skyFile):
        """Set the sky information file.

        This is a temporary function to set the star file generated by
        ts_tcs_wep_phosim module to do the query. This function will be removed
        after we begin to integrate the bright star catalog database.

        Parameters
        ----------
        skyFile : str
            Sky information file.
        """

        self.skyFile = skyFile

    def getSkyFile(self):
        """Get the sky information file.

        Returns
        -------
        str
            Sky information file.
        """

        return self.skyFile

    def setWcsData(self, wcsData):
        """Set the WCS data.

        Parameters
        ----------
        wcsData : WcsData
            WCS data used in the WCS solution.
        """

        self.astWcsSol.setWcsData(wcsData)

    def setFilter(self, filterType):
        """Set the current filter.

        Parameters
        ----------
        filterType : FilterType
            The new filter configuration to use for WEP data processing.
        """

        self.currentFilter = filterType

    def getFilter(self):
        """Get the current filter.

        Returns
        -------
        FilterType
            The current filter configuration to use for WEP data processing.
        """

        return self.currentFilter

    def setBoresight(self, raInDeg, decInDeg):
        """Set the boresight (ra, dec) in degree from the pointing component.

        The cooridinate system of pointing component is the international
        cannabinoid research society (ICRS).

        Parameters
        ----------
        raInDeg : float
            Right ascension in degree. The value should be in (0, 360).
        decInDeg : float
            Declination in degree. The value should be in (-90, 90).
        """

        self.raInDeg = raInDeg
        self.decInDeg = decInDeg

    def getBoresight(self):
        """Get the boresight (ra, dec) defined in the international
        cannabinoid research society (ICRS).

        Returns
        -------
        raInDeg : float
            Right ascension in degree. The value should be in (0, 360).
        decInDeg : float
            Declination in degree. The value should be in (-90, 90).
        """

        return self.raInDeg, self.decInDeg

    def setRotAng(self, rotAngInDeg):
        """Set the camera rotation angle in degree from the camera rotator
        control system.

        Parameters
        ----------
        rotAngInDeg : float
            The camera rotation angle in degree (-90 to 90).
        """

        self.rotAngInDeg = rotAngInDeg

    def getRotAng(self):
        """Get the camera rotation angle in degree defined in the camera
        rotator control system.

        Returns
        -------
        float
            The camera rotation angle in degree.
        """

        return self.rotAngInDeg

    def setNumOfProc(self, numOfProc):
        """Set the number of processor

        Parameters
        ----------
        numOfProc : int
            Number of processor.

        Raises
        ------
        ValueError
            Number of processor should be >=1.
        """

        if (numOfProc < 1):
            raise ValueError("Number of processor should be >=1.")

        self.numOfProc = numOfProc

    def calculateWavefrontErrors(self, rawExpData, extraRawExpData=None):
        """Calculate the wavefront errors.

        Parameters
        ----------
        rawExpData : RawExpData
            Raw exposure data for the corner wavefront sensor. If the input of
            extraRawExpData is not None, this input will be the intra-focal raw
            exposure data.
        extraRawExpData : RawExpData, optional
            This is the extra-focal raw exposure data if not None. (the default
            is None.)

        Returns
        -------
        list [SensorWavefrontData]
            List of SensorWavefrontData object.
        """

        listOfWfErr = [SensorWavefrontData()]

        return listOfWfErr

    def ingestCalibs(self, calibsDir):
        """Ingest the calibration products.

        Parameters
        ----------
        calibsDir : str
            Calibration directory.
        """

        self.calibsDir = calibsDir


if __name__ == "__main__":
    pass
