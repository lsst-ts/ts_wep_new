import unittest

from lsst.ts.wep.Utility import CamType
from lsst.ts.wep.ctrlIntf.WEPCalculationFactory import WEPCalculationFactory
from lsst.ts.wep.ctrlIntf.WEPCalculationOfLsstCam import WEPCalculationOfLsstCam
from lsst.ts.wep.ctrlIntf.WEPCalculationOfLsstFamCam import WEPCalculationOfLsstFamCam
from lsst.ts.wep.ctrlIntf.WEPCalculationOfComCam import WEPCalculationOfComCam


class TestWEPCalculationFactory(unittest.TestCase):
    """Test the WEPCalculationFactory class."""

    def testGetCalculatorOfLsstCam(self):

        calculator = WEPCalculationFactory.getCalculator(CamType.LsstCam)

        self.assertTrue(isinstance(calculator, WEPCalculationOfLsstCam))

    def testGetCalculatorOfLsstFamCam(self):

        calculator = WEPCalculationFactory.getCalculator(CamType.LsstFamCam)

        self.assertTrue(isinstance(calculator, WEPCalculationOfLsstFamCam))

    def testGetCalculatorOfComCam(self):

        calculator = WEPCalculationFactory.getCalculator(CamType.ComCam)

        self.assertTrue(isinstance(calculator, WEPCalculationOfComCam))

    def testGetCalculatorOfWrongCamType(self):

        self.assertRaises(ValueError, WEPCalculationFactory.getCalculator,
                          "wrongInst")


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
