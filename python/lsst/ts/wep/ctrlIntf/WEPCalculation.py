import os

from lsst.ts.wep.Utility import getModulePath, getConfigDir, BscDbType, \
    FilterType
from lsst.ts.wep.CamDataCollector import CamDataCollector
from lsst.ts.wep.CamIsrWrapper import CamIsrWrapper
from lsst.ts.wep.SourceProcessor import SourceProcessor
from lsst.ts.wep.SourceSelector import SourceSelector
from lsst.ts.wep.WfEstimator import WfEstimator
from lsst.ts.wep.WepController import WepController
from lsst.ts.wep.ctrlIntf.SensorWavefrontData import SensorWavefrontData
from lsst.ts.wep.ParamReader import ParamReader


class WEPCalculation(object):
    """Base class for converting the wavefront images into wavefront errors.

    There will be different implementations of this for different
    types of CCDs (normal, full array mode, comcam, cmos, shwfs).
    """

    def __init__(self, astWcsSol, camType, isrDir,
                 settingFileName="default.yaml"):
        """Construct an WEP calculation object.

        Parameters
        ----------
        astWcsSol : AstWcsSol
            AST world coordinate system (WCS) solution.
        camType : enum 'CamType'
            Camera type.
        isrDir : str
            Instrument signature remocal (ISR) directory. This directory will
            have the input and output that the data butler needs.
        settingFileName : str, optional
            Setting file name. (the default is "default.yaml".)
        """

        super().__init__()

        # This attribute is just a stakeholder here since there is no detail of
        # AST WCS solution yet
        self.astWcsSol = astWcsSol

        # ISR directory that the data butler uses
        self.isrDir = isrDir

        # Number of processors for WEP to use
        # This is just a stakeholder at this moment
        self.numOfProc = 1

        # Boresight infomation
        self.raInDeg = 0.0
        self.decInDeg = 0.0

        # Sky rotation angle
        self.rotSkyPos = 0.0

        # Sky information file fro the temporary use
        self.skyFile = ""

        # Default setting file
        settingFilePath = os.path.join(getConfigDir(), settingFileName)
        self.settingFile = ParamReader(settingFilePath)

        # Configure the WEP controller
        self.wepCntlr = self._configWepController(camType, isrDir,
                                                  settingFileName)

    def _configWepController(self, camType, isrDir, settingFileName):
        """Configure the WEP controller.

        WEP: wavefront estimation pipeline.

        Parameters
        ----------
        camType : enum 'CamType'
            Camera type.
        isrDir : str
            Instrument signature remocal (ISR) directory. This directory will
            have the input and output that the data butler needs.
        settingFileName : str
            Setting file name.

        Returns
        -------
        WepController
            Configured WEP controller.
        """

        dataCollector = CamDataCollector(isrDir)
        isrWrapper = CamIsrWrapper(isrDir)

        bscDbType = self._getBscDbType()
        sourSelc = self._configSourceSelector(camType, bscDbType,
                                              settingFileName)

        sourProc = SourceProcessor(settingFileName=settingFileName)
        wfsEsti = self._configWfEstimator(camType)

        wepCntlr = WepController(dataCollector, isrWrapper, sourSelc,
                                 sourProc, wfsEsti)

        return wepCntlr

    def _getBscDbType(self):
        """Get the bright star catalog (BSC) database type.

        Returns
        -------
        enum 'BscDbType'
            BSC database type.

        Raises
        ------
        ValueError
            The bscDb is not supported.
        """

        bscDb = self.settingFile.getSetting("bscDb")
        if (bscDb == "localDb"):
            return BscDbType.LocalDb
        elif (bscDb == "file"):
            return BscDbType.LocalDbForStarFile
        else:
            raise ValueError("The bscDb (%s) is not supported." % bscDb)

    def _configSourceSelector(self, camType, bscDbType, settingFileName):
        """Configue the source selector.

        Parameters
        ----------
        camType : enum 'CamType'
            Camera type.
        bscDbType : enum 'BscDbType'
            Bright star catalog (BSC) database type.
        settingFileName : str
            Setting file name.

        Returns
        -------
        SourceSelector
            Configured source selector.

        Raises
        ------
        ValueError
            WEPCalculation does not support this bscDbType yet.
        """

        sourSelc = SourceSelector(camType, bscDbType,
                                  settingFileName=settingFileName)
        sourSelc.setFilter(FilterType.REF)

        if (bscDbType == BscDbType.LocalDbForStarFile):
            dbAdress = os.path.join(getModulePath(), "tests", "testData",
                                    "bsc.db3")
            sourSelc.connect(dbAdress)
        else:
            raise ValueError("WEPCalculation does not support %s yet." % bscDbType)

        return sourSelc

    def _configWfEstimator(self, camType):
        """Configure the wavefront estimator.

        Returns
        -------
        WfEstimator
            Configured wavefront estimator.
        """

        configDir = getConfigDir()
        instDir = os.path.join(configDir, "cwfs", "instData")
        algoDir = os.path.join(configDir, "cwfs", "algo")
        wfsEsti = WfEstimator(instDir, algoDir)

        solver = self.settingFile.getSetting("poissonSolver")
        opticalModel = self.settingFile.getSetting("opticalModel")
        defocalDisInMm = self.settingFile.getSetting("dofocalDistInMm")
        donutImgSizeInPixel = self.settingFile.getSetting("donutImgSizeInPixel")
        wfsEsti.config(solver=solver, camType=camType,
                       opticalModel=opticalModel,
                       defocalDisInMm=defocalDisInMm,
                       sizeInPix=donutImgSizeInPixel)

        return wfsEsti

    def getWepCntlr(self):
        """Get the configured WEP controller.

        Returns
        -------
        WepController
            Configured WEP controller.
        """

        return self.wepCntlr

    def disconnect(self):
        """Disconnect the database."""

        sourSelc = self.wepCntlr.getSourSelc()
        sourSelc.disconnect()

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
        filterType : enum 'FilterType'
            The new filter configuration to use for WEP data processing.
        """

        sourSelc = self.wepCntlr.getSourSelc()
        sourSelc.setFilter(filterType)

    def getFilter(self):
        """Get the current filter.

        Returns
        -------
        enum 'FilterType'
            The current filter configuration to use for WEP data processing.
        """

        sourSelc = self.wepCntlr.getSourSelc()
        return sourSelc.getFilter()

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

        # In the WCS solution provided by SIMS team, the input angle is sky
        # rotation angle. We do not know its relationship with the camera
        # rotation angle yet.

        self.rotSkyPos = rotAngInDeg

    def getRotAng(self):
        """Get the camera rotation angle in degree defined in the camera
        rotator control system.

        Returns
        -------
        float
            The camera rotation angle in degree.
        """

        # In the WCS solution provided by SIMS team, the input angle is sky
        # rotation angle. We do not know its relationship with the camera
        # rotation angle yet.

        return self.rotSkyPos

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

        # Discuss with Chris that we put this into the configuration file or
        # not

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

        self._genCamMapperIfNeed()

        dataCollector = self.wepCntlr.getDataCollector()

        calibFiles = os.path.join(calibsDir, "*")
        dataCollector.ingestCalibs(calibFiles)

    def _genCamMapperIfNeed(self):
        """Generate the camera mapper file if it is needed.

        The mapper file is used by the data butler.

        Raises
        ------
        ValueError
            Mapper is not supported yet.
        """

        mapperFile = os.path.join(self.isrDir, "_mapper")
        if (not os.path.exists(mapperFile)):

            dataCollector = self.wepCntlr.getDataCollector()

            camMapper = self.settingFile.getSetting("camMapper")
            if (camMapper == "phosim"):
                dataCollector.genPhoSimMapper()
            else:
                raise ValueError("Mapper (%s) is not supported yet." % camMapper)


if __name__ == "__main__":
    pass
