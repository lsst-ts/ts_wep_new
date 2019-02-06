from lsst.ts.wep.ctrlIntf.WEPCalculation import WEPCalculation 
from lsst.ts.wep.ctrlIntf.AstWcsSol import AstWcsSol


class WEPCalculationOfLsstCam(WEPCalculation):
    """The concrete child class of WEPCalculation of the LSST camera (corner
    wavefront sensor)."""

    def __init__(self):
        """Construct an WEP calculation of LSST camera object."""
        super(WEPCalculationOfLsstCam, self).__init__(AstWcsSol())

        self.visit = 0
        self.snap = 0
        self.rawExpDir = ""

    def ingestRawExp(self, visit, snap, rawExpDir):
        """Ingest the raw exposure data. The raw exposures should be in the
        local directory already.

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

        self.visit = visit
        self.snap = snap
        self.rawExpDir = rawExpDir

    def _resetRawExpInfo(self):
        """Reset the raw exposure information."""

        self.visit = 0
        self.snap = 0
        self.rawExpDir = ""


if __name__ == "__main__":
    pass
