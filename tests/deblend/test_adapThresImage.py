import os
import numpy as np
import unittest

from lsst.ts.wep.deblend.AdapThresImage import AdapThresImage
from lsst.ts.wep.Utility import getModulePath


class TestAdapThresImage(unittest.TestCase):
    """Test the AdapThresImage class."""

    def setUp(self):

        # Get the path of module
        modulePath = getModulePath()

        # Image file
        imageFolderPath = os.path.join(modulePath, "tests", "testData",
                                       "testImages", "LSST_NE_SN25")
        imageName = 'z11_0.25_intra.txt'
        imageFile = os.path.join(imageFolderPath, imageName)

        # Set the image file
        self.adapImage = AdapThresImage()
        self.adapImage.setImg(imageFile=imageFile)

    def testFunc(self):

        # Answer of centroid
        ansCx = 61
        ansCy = 61
        ansR = 38

        # Allowed tolerace
        tor = 5

        # Calculate the centroid
        realCx, realCy, realR, imgBinary = self.adapImage.getCenterAndR_adap()

        delta = self._calcSqrtDelta(realCx, realCy, realR, ansCx, ansCy, ansR)
        self.assertLess(delta, tor)

        realCx, realCy, realR = self.adapImage.getCenterAndR()
        delta = self._calcSqrtDelta(realCx, realCy, realR, ansCx, ansCy, ansR)
        self.assertLess(delta, tor)

    def _calcSqrtDelta(self, realCx, realCy, realR, ansCx, ansCy, ansR):

        delta = np.sqrt((realCx-ansCx)**2 + (realCy-ansCy)**2 + (realR-ansR)**2)

        return delta


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
