import os
import numpy as np
import unittest

from lsst.ts.wep.cwfs.Image import Image
from lsst.ts.wep.Utility import getModulePath


class TestImage(unittest.TestCase):
    """Test the Image class."""

    def setUp(self):

        # Get the path of module
        modulePath = getModulePath()

        # Define the image folder and image names
        # Image data -- Don't know the final image format.
        # It is noted that image.readFile inuts is based on the txt file
        imageFolderPath = os.path.join(modulePath, "tests", "testData",
                                       "testImages", "LSST_NE_SN25")
        imgName = "z11_0.25_intra.txt"

        # Image files Path
        imgFile = os.path.join(imageFolderPath, imgName)

        # There is the difference between intra and extra images
        self.img = Image()
        self.img.setImg(imageFile=imgFile)

    def testZeroImg(self):

        # Creat a zero image
        zeroImg = Image()
        zeroImg.setImg(image=np.zeros([4, 4]))
        self.assertEqual(np.sum(zeroImg.image), 0)

        # Update Image
        zeroImg.updateImage(np.ones([4, 4]))
        self.assertEqual(np.sum(zeroImg.image), 16)

        realcx, realcy, realR, imgBinary = \
            zeroImg.getCenterAndR_ef(randNumFilePath=None, checkEntropy=True)

        self.assertEqual(realcx, [])

        # update to the random image
        zeroImg.updateImage(np.random.rand(100, 100))
        realcx, realcy, realR, imgBinary = \
            zeroImg.getCenterAndR_ef(randNumFilePath=None, checkEntropy=True)
        self.assertEqual(realcx, [])

    def testImg(self):

        realcx, realcy, realR, imgBinary = \
            self.img.getCenterAndR_ef(randNumFilePath=None, checkEntropy=True)
        self.assertEqual(int(realcx), 61)
        self.assertEqual(int(realcy), 61)
        self.assertGreater(int(realR), 35)

        # Calculate the S/N
        # Add the noise to the image
        noisedImg = self.img.image + np.random.random(self.img.image.shape)*0.1
        self.img.setImg(image=noisedImg)
        snr = self.img.getSNR()

        self.assertGreater(snr, 15)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
