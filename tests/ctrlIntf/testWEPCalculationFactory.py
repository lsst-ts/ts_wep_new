import unittest

from lsst.ts.wep.Utility import CamType
from lsst.ts.wep.ctrlIntf.WEPCalculationFactory import WEPCalculationFactory
from lsst.ts.wep.ctrlIntf.WEPCalculationOfLsstCam import WEPCalculationOfLsstCam
from lsst.ts.wep.ctrlIntf.WEPCalculationOfLsstFamCam import WEPCalculationOfLsstFamCam
from lsst.ts.wep.ctrlIntf.WEPCalculationOfComCam import WEPCalculationOfComCam


class TestWEPCalculationFactory(unittest.TestCase):
    """Test the WEPCalculationFactory class."""

    def setUp(self):

        self.wepCalculationFactory = WEPCalculationFactory()

    def testGetCalculatorOfLsstCam(self):

        calculator = self.wepCalculationFactory.getCalculator(CamType.LsstCam)

        self.assertTrue(isinstance(calculator, WEPCalculationOfLsstCam))

    def testGetCalculatorOfLsstFamCam(self):

        calculator = self.wepCalculationFactory.getCalculator(
            CamType.LsstFamCam)

        self.assertTrue(isinstance(calculator, WEPCalculationOfLsstFamCam))

    def testGetCalculatorOfComCam(self):

        calculator = self.wepCalculationFactory.getCalculator(CamType.ComCam)

        self.assertTrue(isinstance(calculator, WEPCalculationOfComCam))

    def testGetCalculatorOfWrongCamType(self):

        self.assertRaises(ValueError, self.wepCalculationFactory.getCalculator,
                          "wrongInst")


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
