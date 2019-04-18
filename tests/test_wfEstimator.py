import os
import numpy as np
import unittest

from lsst.ts.wep.WfEstimator import WfEstimator
from lsst.ts.wep.Utility import getModulePath, getConfigDir, DefocalType


class TestWfEsitmator(unittest.TestCase):
    """Test the wavefront estimator class."""

    def setUp(self):

        cwfsConfigDir = os.path.join(getConfigDir(), "cwfs")
        instDir = os.path.join(cwfsConfigDir, "instData")
        algoDir = os.path.join(cwfsConfigDir, "algo")
        self.wfsEst = WfEstimator(instDir, algoDir)

        # Define the image folder and image names
        # It is noted that image.readFile inuts is based on the txt file.
        self.modulePath = getModulePath()
        imageFolderPath = os.path.join(self.modulePath, "tests", "testData",
                                       "testImages", "LSST_NE_SN25")
        intra_image_name = "z11_0.25_intra.txt"
        extra_image_name = "z11_0.25_extra.txt"

        # Path to image files
        self.intraImgFile = os.path.join(imageFolderPath, intra_image_name)
        self.extraImgFile = os.path.join(imageFolderPath, extra_image_name)

        # Field XY position
        self.fieldXY = (1.185, 1.185)

    def testCalWfsErrOfExp(self):

        # Setup the images
        self.wfsEst.setImg(self.fieldXY, DefocalType.Intra,
                           imageFile=self.intraImgFile)
        self.wfsEst.setImg(self.fieldXY, DefocalType.Extra,
                           imageFile=self.extraImgFile)

        # Test the images are set.
        self.assertEqual(self.wfsEst.ImgIntra.defocalType, DefocalType.Intra)
        self.assertEqual(self.wfsEst.ImgExtra.defocalType, DefocalType.Extra)

        # Setup the configuration

        # Try to catch the error
        try:
            self.wfsEst.config(solver="exp", defocalDisInMm=3, debugLevel=0)
        except ValueError:
            print("Catch the wrong instrument.")

        # If the configuration is reset, the images are needed to be set again.
        self.wfsEst.config(solver="exp", instName="lsst",
                           opticalModel="offAxis", defocalDisInMm=1.0,
                           sizeInPix=120, debugLevel=0)

        # Evaluate the wavefront error
        wfsError = [2.593, 14.102, -8.470, 3.676, 1.467, -9.724, 8.207,
                    -192.839, 0.978, 1.568, 4.197, -0.391, 1.551, 1.235,
                    -1.699, 2.140, -0.296, -2.113, 1.188]
        zer4UpNm = self.wfsEst.calWfsErr()
        self.assertAlmostEqual(np.sum(np.abs(zer4UpNm-np.array(wfsError))), 0,
                               places=1)

    def testCalWfsErrOfFft(self):

        # Reset the wavefront images
        self.wfsEst.setImg(self.fieldXY, DefocalType.Intra,
                           imageFile=self.intraImgFile)
        self.wfsEst.setImg(self.fieldXY, DefocalType.Extra,
                           imageFile=self.extraImgFile)

        # Change the algorithm to fft
        self.wfsEst.config(solver="fft", instName="lsst",
                           opticalModel="offAxis", defocalDisInMm=1.0,
                           sizeInPix=120, debugLevel=0)

        # Evaluate the wavefront error
        wfsError = [12.484, 10.358, -6.674, -0.043, -1.768, -15.593, 12.511,
                    -192.382, 0.195, 4.074, 9.577, -1.930, 3.538, 3.420,
                    -3.610, 3.547, -0.679, -2.943, 1.101]
        zer4UpNm = self.wfsEst.calWfsErr()
        self.assertAlmostEqual(np.sum(np.abs(zer4UpNm-np.array(wfsError))), 0,
                               places=1)

        # Test to reset the data
        self.wfsEst.reset()
        self.assertEqual(np.sum(self.wfsEst.getAlgo().getZer4UpInNm()), 0)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
