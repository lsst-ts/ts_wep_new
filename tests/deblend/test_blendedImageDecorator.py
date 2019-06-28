import os
import numpy as np
import unittest

from lsst.ts.wep.deblend.AdapThresImage import AdapThresImage
from lsst.ts.wep.deblend.BlendedImageDecorator import BlendedImageDecorator
from lsst.ts.wep.Utility import getModulePath


class TestBlendedImageDecorator(unittest.TestCase):
    """Test the BlendedImageDecorator class."""

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

    def testBlendedImage(self):

        # Generate the blended image
        image, imageMain, imageNeighbor, neighborX, neighborY = \
            self.adapImage.generateMultiDonut(1.3, 0.1, 0.0)

        blendImage = BlendedImageDecorator()
        blendImage.setImg(image=image)

        # Do the deblending
        imgDeblend, realcx, realcy = \
            blendImage.deblendDonut((neighborX, neighborY))

        # Do the comparison
        delta = np.sum(np.abs(imageMain-imgDeblend))
        diffRatio = delta/np.sum(np.abs(imageMain))
        self.assertLess(diffRatio, 0.01)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
