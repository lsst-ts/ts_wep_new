from lsst.ts.wep.ctrlIntf.WEPCalculation import WEPCalculation


class WEPCalculationOfPiston(WEPCalculation):
    """The child class of WEPCalculation that gets the defocal images by the
    camera piston."""

    def __init__(self, astWcsSol):
        """Construct an WEP calculation of piston object."""
        super(WEPCalculationOfPiston, self).__init__(astWcsSol)

        self.intraVisit = []
        self.extraVisit = []

        self.intraSnap = []
        self.extraSnap = []

        self.intraRawExpDir = []
        self.extraRawExpDir = []

        self.defocalDisInMm = 1.5

    def ingestIntraRawExp(self, visit, snap, rawExpDir):
        """Ingest the intra-focal raw exposure data. The raw exposures should
        be in the local directory already.

        Parameters
        ----------
        visit : int
            Unique visit Id.
        snap : int
            Snap (0 or 1).
        rawExpDir : str
            Raw exposure directory in the local disk.
        """

        # Ingest the images

        self.intraVisit.append(visit)
        self.intraSnap.append(snap)
        self.intraRawExpDir.append(rawExpDir)

    def ingestExtraRawExp(self, visit, snap, rawExpDir):
        """Ingest the extra-focal raw exposure data. The raw exposures should
        be in the local directory already.

        Parameters
        ----------
        visit : int
            Unique visit Id.
        snap : int
            Snap (0 or 1).
        rawExpDir : str
            Raw exposure directory in the local disk.
        """

        # Ingest the images

        self.extraVisit.append(visit)
        self.extraSnap.append(snap)
        self.extraRawExpDir.append(rawExpDir)

    def _resetRawExpInfo(self):
        """Reset the raw exposure information."""

        self.intraVisit = []
        self.extraVisit = []

        self.intraSnap = []
        self.extraSnap = []

        self.intraRawExpDir = []
        self.extraRawExpDir = []

    def setDefocalDisInMm(self, defocalDisInMm):
        """Set the defocal distance in mm.

        Parameters
        ----------
        float
            Defocal distance in mm.
        """

        self.defocalDisInMm = defocalDisInMm

    def getDefocalDisInMm(self):
        """Set the defocal distance in mm.

        Returns
        -------
        float
            Defocal distance in mm.
        """

        return self.defocalDisInMm


if __name__ == "__main__":
    pass
