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
        self.assertEqual(defocalDisInMm, 1.5225)

    def testSetDefocalDisInMm(self):

        defocalDisInMm = 2.0
        self.wepCalculationOfPiston.setDefocalDisInMm(defocalDisInMm)

        self.assertAlmostEqual(self.wepCalculationOfPiston.getDefocalDisInMm(),
                               2.03)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
