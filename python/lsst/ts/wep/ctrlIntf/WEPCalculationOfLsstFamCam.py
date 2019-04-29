from lsst.ts.wep.ctrlIntf.WEPCalculationOfPiston import WEPCalculationOfPiston
from lsst.ts.wep.ctrlIntf.AstWcsSol import AstWcsSol


class WEPCalculationOfLsstFamCam(WEPCalculationOfPiston):
    """The concrete child class of WEPCalculationOfPiston of the LSST
    full-array mode (FAM) camera."""

    def __init__(self, isrDir):
        """Construct an WEP calculation of LSST FAM camera object.

        Parameters
        ----------
        isrDir : str
            Instrument signature remocal (ISR) directory. This directory will
            have the input and output that the data butler needs.
        """
        super(WEPCalculationOfLsstFamCam, self).__init__(AstWcsSol(), isrDir)


if __name__ == "__main__":
    pass
