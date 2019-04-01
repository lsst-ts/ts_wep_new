from lsst.ts.wep.ctrlIntf.WEPCalculationOfPiston import WEPCalculationOfPiston
from lsst.ts.wep.ctrlIntf.AstWcsSol import AstWcsSol


class WEPCalculationOfComCam(WEPCalculationOfPiston):
    """The concrete child class of WEPCalculationOfPiston of the commionning
    camera (ComCam)."""

    def __init__(self):
        """Construct an WEP calculation of ComCam object."""
        super(WEPCalculationOfComCam, self).__init__(AstWcsSol())


if __name__ == "__main__":
    pass
