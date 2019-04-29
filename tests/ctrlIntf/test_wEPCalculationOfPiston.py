import unittest

from lsst.ts.wep.Utility import CamType
from lsst.ts.wep.ctrlIntf.WEPCalculationOfPiston import WEPCalculationOfPiston
from lsst.ts.wep.ctrlIntf.AstWcsSol import AstWcsSol


class TestWEPCalculationOfPiston(unittest.TestCase):
    """Test the WEPCalculationOfPiston class."""

    def setUp(self):

        self.wepCalculationOfPiston = WEPCalculationOfPiston(
            AstWcsSol(), CamType.ComCam, "")

    def testGetDefocalDisInMm(self):

        defocalDisInMm = self.wepCalculationOfPiston.getDefocalDisInMm()
        self.assertEqual(defocalDisInMm,
                         WEPCalculationOfPiston.DEFOCAL_DIS_IN_MM)

    def testSetDefocalDisInMm(self):

        defocalDisInMm = 3.0
        self.wepCalculationOfPiston.setDefocalDisInMm(defocalDisInMm)

        self.assertEqual(self.wepCalculationOfPiston.getDefocalDisInMm(),
                         defocalDisInMm)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
