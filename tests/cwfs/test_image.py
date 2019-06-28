import os
import numpy as np
import unittest

from lsst.ts.wep.cwfs.Image import Image
from lsst.ts.wep.cwfs.CentroidRandomWalk import CentroidRandomWalk
from lsst.ts.wep.Utility import getModulePath


class TestImage(unittest.TestCase):
    """Test the Image class."""

    def setUp(self):

        self.testDataDir = os.path.join(getModulePath(), "tests", "testData")
        self.imgFile = os.path.join(self.testDataDir, "testImages",
                                    "LSST_NE_SN25", "z11_0.25_intra.txt")

        self.img = Image()
        self.img.setImg(imageFile=self.imgFile)

    def testGetCentroidFind(self):

        centroidFind = self.img.getCentroidFind()
        self.assertTrue(isinstance(centroidFind, CentroidRandomWalk))

    def testGetImg(self):

        img = self.img.getImg()
        self.assertEqual(img.shape, (120, 120))

    def testGetImgFilePath(self):

        imgFilePath = self.img.getImgFilePath()
        self.assertEqual(imgFilePath, self.imgFile)

    def testSetImgByImageArray(self):

        newImg = np.random.rand(5, 5)
        self.img.setImg(image=newImg)

        self.assertTrue(np.all(self.img.getImg() == newImg))
        self.assertEqual(self.img.getImgFilePath(), "")

    def testSetImgByFitsFile(self):

        opdFitsFile = os.path.join(self.testDataDir, "opdOutput", "9005000",
                                   "opd_9005000_0.fits.gz")
        self.img.setImg(imageFile=opdFitsFile)

        img = self.img.getImg()
        self.assertEqual(img.shape, (255, 255))

        imgFilePath = self.img.getImgFilePath()
        self.assertEqual(imgFilePath, opdFitsFile)

    def testUpdateImage(self):

        newImg = np.random.rand(5, 5)
        self.img.updateImage(newImg)

        self.assertTrue(np.all(self.img.getImg() == newImg))

    def testUpdateImageWithNoHoldImage(self):

        img = Image()

        newImg = np.random.rand(5, 5)
        self.assertWarns(UserWarning, img.updateImage, newImg)

    def testGetCenterAndR_ef(self):

        realcx, realcy, realR = self.img.getCenterAndR()
        self.assertEqual(int(realcx), 61)
        self.assertEqual(int(realcy), 61)
        self.assertGreater(int(realR), 35)

    def testGetSNR(self):

        # Add the noise to the image
        image = self.img.getImg()
        noisedImg = image + np.random.random(image.shape) * 0.1

        self.img.setImg(image=noisedImg)
        snr = self.img.getSNR()

        self.assertGreater(snr, 15)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
