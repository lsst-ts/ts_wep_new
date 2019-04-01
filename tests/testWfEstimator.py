import os
import numpy as np
import unittest

from lsst.ts.wep.WfEstimator import WfEstimator
from lsst.ts.wep.Utility import getModulePath


class TestWfEsitmator(unittest.TestCase):
    """Test the wavefront estimator class."""

    def setUp(self):

        # Get the path of module
        self.modulePath = getModulePath()

        # Define the instrument folder
        instruFolderPath = os.path.join(self.modulePath, "configData", "cwfs",
                                        "instruData")

        # Define the algorithm folder
        algoFolderPath = os.path.join(self.modulePath, "configData", "cwfs",
                                      "algo")

        # Decalre the WfEstimator
        self.wfsEst = WfEstimator(instruFolderPath, algoFolderPath)

    def testFunc(self):

        # Define the image folder and image names
        # Image data -- Don't know the final image format.
        # It is noted that image.readFile inuts is based on the txt file.
        imageFolderPath = os.path.join(self.modulePath, "tests", "testData",
                                       "testImages", "LSST_NE_SN25")
        intra_image_name = "z11_0.25_intra.txt"
        extra_image_name = "z11_0.25_extra.txt"

        # Path to image files
        intraImgFile = os.path.join(imageFolderPath, intra_image_name)
        extraImgFile = os.path.join(imageFolderPath, extra_image_name)

        # Field XY position
        fieldXY = [1.185, 1.185]

        # Setup the images
        self.wfsEst.setImg(fieldXY, imageFile=intraImgFile,
                           defocalType="intra")
        self.wfsEst.setImg(fieldXY, imageFile=extraImgFile,
                           defocalType="extra")

        # Test the images are set.
        self.assertEqual(self.wfsEst.ImgIntra.atype,
                         self.wfsEst.ImgIntra.INTRA)
        self.assertEqual(self.wfsEst.ImgExtra.atype,
                         self.wfsEst.ImgExtra.EXTRA)

        # Setup the configuration

        # Try to catch the error
        try:
            self.wfsEst.config(solver="exp", defocalDisInMm=3, debugLevel=0)
        except ValueError:
            print("Catch the wrong instrument.")

        # If the configuration is reset, the images are needed to be set again.
        self.wfsEst.config(solver="exp", instName="lsst", sizeInPix=120,
                           opticalModel="offAxis", debugLevel=0)

        # Test the setting of algorithm and instrument
        self.assertEqual(self.wfsEst.inst.instName, "lsst")
        self.assertEqual(self.wfsEst.algo.algoName, "exp")

        # Evaluate the wavefront error
        wfsError = [2.593, 14.102, -8.470, 3.676, 1.467, -9.724, 8.207,
                    -192.839, 0.978, 1.568, 4.197, -0.391, 1.551, 1.235,
                    -1.699, 2.140, -0.296, -2.113, 1.188]
        zer4UpNm = self.wfsEst.calWfsErr()
        self.assertAlmostEqual(np.sum(np.abs(zer4UpNm-np.array(wfsError))), 0,
                               places=1)

        # Reset the wavefront images
        self.wfsEst.setImg(fieldXY, imageFile=intraImgFile,
                           defocalType="intra")
        self.wfsEst.setImg(fieldXY, imageFile=extraImgFile,
                           defocalType="extra")

        # Change the algorithm to fft
        self.wfsEst.config(solver="fft")
        self.assertEqual(self.wfsEst.algo.algoName, "fft")

        # Evaluate the wavefront error
        wfsError = [12.484, 10.358, -6.674, -0.043, -1.768, -15.593, 12.511,
                    -192.382, 0.195, 4.074, 9.577, -1.930, 3.538, 3.420,
                    -3.610, 3.547, -0.679, -2.943, 1.101]
        zer4UpNm = self.wfsEst.calWfsErr()
        self.assertAlmostEqual(np.sum(np.abs(zer4UpNm-np.array(wfsError))), 0,
                               places=1)

        # Test to output the parameters
        filename = "outputParameter"
        self.wfsEst.outParam(filename=filename)
        self.assertTrue(os.path.isfile(filename))
        os.remove(filename)

        # Test to reset the data
        self.wfsEst.reset()
        self.assertEqual(np.sum(self.wfsEst.algo.zer4UpNm), 0)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
