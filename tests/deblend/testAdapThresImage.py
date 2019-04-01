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

        self.zeroImage = AdapThresImage()
        self.zeroImage.setImg(image=np.zeros([120, 120]))

        self.randImage = AdapThresImage()
        self.randImage.setImg(image=np.random.rand(120, 120))

    def testFunc(self):

        # Calculate the centroid
        realcx, realcy, realR, imgBinary = self.adapImage.getCenterAndR_adap()
        self.assertEqual([round(realcx), round(realcy), round(realR)],
                         [61, 61, 38])

        realcx, realcy, realR, imgBinary = \
            self.adapImage.getCenterAndR_ef(checkEntropy=True)
        self.assertEqual([round(realcx), round(realcy), round(realR)],
                         [61, 61, 38])

        realcx, realcy, realR, imgBinary = \
            self.zeroImage.getCenterAndR_ef(checkEntropy=True)
        self.assertEqual(realcx, [])

        realcx, realcy, realR, imgBinary = \
            self.randImage.getCenterAndR_ef(checkEntropy=True)
        self.assertEqual(realcx, [])


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
