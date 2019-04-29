from lsst.ts.wep.Utility import FilterType

from lsst.ts.wep.ctrlIntf.SensorWavefrontData import SensorWavefrontData


class WEPCalculation(object):
    """Base class for converting the wavefront images into wavefront errors.

    There will be different implementations of this for different
    types of CCDs (normal, full array mode, comcam, cmos, shwfs).
    """

    def __init__(self, astWcsSol, isrDir):
        """Construct an WEP calculation object.

        Parameters
        ----------
        astWcsSol : AstWcsSol
            AST world coordinate system (WCS) solution.
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

        self.currentFilter = FilterType.REF
        self.raInDeg = 0.0
        self.decInDeg = 0.0
        self.rotAngInDeg = 0.0

        self.calibsDir = ""

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
