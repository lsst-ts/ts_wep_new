from lsst.ts.wep.ctrlIntf.WEPCalculation import WEPCalculation
from lsst.ts.wep.ctrlIntf.AstWcsSol import AstWcsSol


class WEPCalculationOfLsstCam(WEPCalculation):
    """The concrete child class of WEPCalculation of the LSST camera (corner
    wavefront sensor)."""

    def __init__(self):
        """Construct an WEP calculation of LSST camera object."""
        super(WEPCalculationOfLsstCam, self).__init__(AstWcsSol())


if __name__ == "__main__":
    pass
