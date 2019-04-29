from lsst.ts.wep.Utility import CamType
from lsst.ts.wep.ctrlIntf.WEPCalculationOfPiston import WEPCalculationOfPiston
from lsst.ts.wep.ctrlIntf.AstWcsSol import AstWcsSol


class WEPCalculationOfComCam(WEPCalculationOfPiston):
    """The concrete child class of WEPCalculationOfPiston of the commionning
    camera (ComCam)."""

    def __init__(self, isrDir):
        """Construct an WEP calculation of ComCam object.

        Parameters
        ----------
        isrDir : str
            Instrument signature remocal (ISR) directory. This directory will
            have the input and output that the data butler needs.
        """
        super(WEPCalculationOfComCam, self).__init__(AstWcsSol(),
                                                     CamType.ComCam, isrDir)


if __name__ == "__main__":
    pass
