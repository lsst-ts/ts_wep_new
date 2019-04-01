from lsst.ts.wep.ctrlIntf.WEPCalculationOfPiston import WEPCalculationOfPiston
from lsst.ts.wep.ctrlIntf.AstWcsSol import AstWcsSol


class WEPCalculationOfLsstFamCam(WEPCalculationOfPiston):
    """The concrete child class of WEPCalculationOfPiston of the LSST
    full-array mode (FAM) camera."""

    def __init__(self):
        """Construct an WEP calculation of LSST FAM camera object."""
        super(WEPCalculationOfLsstFamCam, self).__init__(AstWcsSol())


if __name__ == "__main__":
    pass
