import os
import numpy as np
import unittest

from lsst.ts.wep.DonutImageCheck import DonutImageCheck
from lsst.ts.wep.Utility import getModulePath


class TestDonutImageCheck(unittest.TestCase):
    """Test the DonutImageCheck class."""

    def setUp(self):

        self.donutImgCheck = DonutImageCheck()

    def testIsEffDonutWithEffImg(self):

        imgFile = os.path.join(getModulePath(), "tests", "testData",
                               "testImages", "LSST_NE_SN25",
                               "z11_0.25_intra.txt")
        donutImg = np.loadtxt(imgFile)

        self.assertTrue(self.donutImgCheck.isEffDonut(donutImg))

    def testIsEffDonutWithZeroImg(self):

        donutImg = np.zeros((120, 120))
        self.assertFalse(self.donutImgCheck.isEffDonut(donutImg))

    def testIsEffDonutWithRandImg(self):

        donutImg = np.random.rand(120, 120)
        self.assertFalse(self.donutImgCheck.isEffDonut(donutImg))


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
