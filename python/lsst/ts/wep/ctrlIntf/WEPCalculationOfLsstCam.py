from lsst.ts.wep.ctrlIntf.WEPCalculation import WEPCalculation
from lsst.ts.wep.ctrlIntf.AstWcsSol import AstWcsSol


class WEPCalculationOfLsstCam(WEPCalculation):
    """The concrete child class of WEPCalculation of the LSST camera (corner
    wavefront sensor)."""

    def __init__(self, isrDir):
        """Construct an WEP calculation of LSST camera object.

        Parameters
        ----------
        isrDir : str
            Instrument signature remocal (ISR) directory. This directory will
            have the input and output that the data butler needs.
        """
        super(WEPCalculationOfLsstCam, self).__init__(AstWcsSol(), isrDir)


if __name__ == "__main__":
    pass
