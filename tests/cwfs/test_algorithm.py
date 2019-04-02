import os
import numpy as np
import unittest

from lsst.ts.wep.cwfs.Instrument import Instrument
from lsst.ts.wep.cwfs.CompensationImageDecorator import CompensationImageDecorator
from lsst.ts.wep.cwfs.Algorithm import Algorithm
from lsst.ts.wep.Utility import getModulePath


class TestAlgorithm(unittest.TestCase):
    """Test the Algorithm class."""

    def setUp(self):

        # Get the path of module
        self.modulePath = getModulePath()

        # Define the instrument folder
        instruFolder = os.path.join(self.modulePath, "configData", "cwfs",
                                    "instruData")

        # Define the algorithm folder
        self.algoFolderPath = os.path.join(self.modulePath, "configData", "cwfs",
                                           "algo")

        # Define the instrument name
        instruName = "lsst"

        # Define the image folder and image names
        # Image data -- Don't know the final image format.
        # It is noted that image.readFile inuts is based on the txt file
        imageFolderPath = os.path.join(self.modulePath, "tests", "testData",
                                       "testImages", "LSST_NE_SN25")
        intra_image_name = "z11_0.25_intra.txt"
        extra_image_name = "z11_0.25_extra.txt"

        # Define fieldXY: [1.185, 1.185] or [0, 0]
        # This is the position of donut on the focal plane in degree
        fieldXY = [1.185, 1.185]

        # Define the optical model: "paraxial", "onAxis", "offAxis"
        self.opticalModel = "offAxis"

        # Image files Path
        intra_image_file = os.path.join(imageFolderPath, intra_image_name)
        extra_image_file = os.path.join(imageFolderPath, extra_image_name)

        # Theree is the difference between intra and extra images
        # I1: intra_focal images, I2: extra_focal Images
        # self.I1 = Image.Image()
        # self.I2 = Image.Image()

        self.I1 = CompensationImageDecorator()
        self.I2 = CompensationImageDecorator()

        self.I1.setImg(fieldXY, imageFile=intra_image_file, atype="intra")
        self.I2.setImg(fieldXY, imageFile=extra_image_file, atype="extra")

        self.inst = Instrument(instruFolder)
        self.inst.config(instruName, self.I1.sizeinPix)

    def testExp(self):

        # Define the algorithm being used: "exp" or "fft"
        useAlgorithm = "exp"

        # Define the algorithm to be used.
        algo = Algorithm(self.algoFolderPath)
        algo.config(useAlgorithm, self.inst, debugLevel=3)
        algo.setDebugLevel(0)
        self.assertEqual(algo.debugLevel, 0)

        # Test functions: itr0() and nextItr()
        algo.itr0(self.inst, self.I1, self.I2, self.opticalModel)
        tmp1 = algo.zer4UpNm
        algo.nextItr(self.inst, self.I1, self.I2, self.opticalModel, nItr=2)
        algo.itr0(self.inst, self.I1, self.I2, self.opticalModel)
        tmp2 = algo.zer4UpNm

        difference = np.sum(np.abs(tmp1-tmp2))
        self.assertEqual(difference, 0)

        # Run it
        algo.runIt(self.inst, self.I1, self.I2, self.opticalModel, tol=1e-3)

        # Check the value
        Zk = algo.zer4UpNm
        self.assertEqual(int(Zk[7]), -192)

        # Reset and check the calculation again
        fieldXY = [self.I1.fieldX, self.I1.fieldY]
        self.I1.setImg(fieldXY, image=self.I1.image0, atype=self.I1.atype)
        self.I2.setImg(fieldXY, image=self.I2.image0, atype=self.I2.atype)
        algo.reset()
        algo.runIt(self.inst, self.I1, self.I2, self.opticalModel, tol=1e-3)
        Zk = algo.zer4UpNm
        self.assertEqual(int(Zk[7]), -192)

    def testFFT(self):

        # Define the algorithm being used: "exp" or "fft"
        useAlgorithm = "fft"

        # Define the algorithm to be used.
        algo = Algorithm(self.algoFolderPath)
        algo.config(useAlgorithm, self.inst, debugLevel=0)

        # Run it
        algo.runIt(self.inst, self.I1, self.I2, self.opticalModel, tol=1e-3)

        # Check the value
        Zk = algo.zer4UpNm
        self.assertEqual(int(Zk[7]), -192)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
