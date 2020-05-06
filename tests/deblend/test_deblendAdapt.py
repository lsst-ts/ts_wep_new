import os
import unittest
import numpy as np

from lsst.ts.wep.Utility import getModulePath
from lsst.ts.wep.deblend.DeblendAdapt import DeblendAdapt


class TestDeblendAdapt(unittest.TestCase):
    """Test the DeblendAdapt class."""

    def setUp(self):

        self.deblend = DeblendAdapt()

    def testDeblendDonutMoreThanOneStarNgr(self):

        iniGuessXY = [(1, 2), (3, 4)]
        self.assertRaises(ValueError, self.deblend.deblendDonut, [], iniGuessXY)

    def testDeblendDonut(self):

        template, imgToDeblend, iniGuessXY = self._genBlendedImg()
        imgDeblend, realcx, realcy = self.deblend.deblendDonut(imgToDeblend,
                                                               iniGuessXY)

        difference = np.sum(np.abs(np.sum(template) - np.sum(imgDeblend)))
        self.assertLess(difference, 20)

        self.assertEqual(np.rint(realcx), 96)
        self.assertEqual(np.rint(realcy), 93)

    def _genBlendedImg(self):

        imageFilePath = os.path.join(getModulePath(), "tests", "testData",
                                     "testImages", "LSST_NE_SN25",
                                     "z11_0.25_intra.txt")
        template = np.loadtxt(imageFilePath)

        image, imageMain, imageNeighbor, neighborX, neighborY = \
            self.deblend.generateMultiDonut(template, 1.3, 0.1, 45.0)

        return template, image, [(neighborX, neighborY)]


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
