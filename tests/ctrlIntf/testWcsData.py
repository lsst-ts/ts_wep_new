import unittest
import numpy as np

from lsst.ts.wep.ctrlIntf.WcsData import WcsData


class TestWcsData(unittest.TestCase):
    """Test the WcsData class."""

    def setUp(self):

        self.wcsCoef = np.arange(4)
        self.wcsData = WcsData(self.wcsCoef)

    def testGetWcsCoef(self):

        wcsCoef = self.wcsData.getWcsCoef()

        delta = np.sum(np.abs(wcsCoef - self.wcsCoef))
        self.assertEqual(delta, 0)

    def testSetWcsCoef(self):

        wcsCoef = np.arange(5)
        self.wcsData.setWcsCoef(wcsCoef)

        wcsCoefInWcsData = self.wcsData.getWcsCoef()
        delta = np.sum(np.abs(wcsCoefInWcsData - wcsCoef))
        self.assertEqual(delta, 0)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
