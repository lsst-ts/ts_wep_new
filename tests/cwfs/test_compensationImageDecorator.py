import os
import numpy as np
import unittest

from lsst.ts.wep.cwfs.Instrument import Instrument
from lsst.ts.wep.cwfs.CompensationImageDecorator import CompensationImageDecorator
from lsst.ts.wep.Utility import getModulePath


class TempAlgo(object):
    """Temporary algorithm class used for the testing."""

    def __init__(self):
        self.parameter = None


class CompensationImageDecoratorTest(unittest.TestCase):
    """Test the CompensationImageDecorator class."""

    def setUp(self):

        # Get the path of module
        modulePath = getModulePath()

        # Define the instrument folder
        instruFolder = os.path.join(modulePath, "configData", "cwfs",
                                    "instruData")

        # Define the instrument name
        instruName = "lsst"
        sensorSamples = 120

        self.inst = Instrument(instruFolder)
        self.inst.config(instruName, sensorSamples)

        # Define the image folder and image names
        # Image data -- Don't know the final image format.
        # It is noted that image.readFile inuts is based on the txt file
        imageFolderPath = os.path.join(modulePath, "tests", "testData",
                                       "testImages", "LSST_NE_SN25")
        intra_image_name = "z11_0.25_intra.txt"
        extra_image_name = "z11_0.25_extra.txt"
        self.imgFilePathIntra = os.path.join(imageFolderPath, intra_image_name)
        self.imgFilePathExtra = os.path.join(imageFolderPath, extra_image_name)

        # Define fieldXY: [1.185, 1.185] or [0, 0]
        # This is the position of donut on the focal plane in degree
        self.fieldXY = [1.185, 1.185]

        # Define the optical model: "paraxial", "onAxis", "offAxis"
        self.opticalModel = "offAxis"

        # Get the true Zk
        zcAnsFilePath = os.path.join(modulePath, "tests", "testData",
                                     "testImages", "validation",
                                     "LSST_NE_SN25_z11_0.25_exp.txt")
        self.zcCol = np.loadtxt(zcAnsFilePath)

    def testFunc(self):

        # Declare the CompensationImageDecorator
        wfsImg = CompensationImageDecorator()

        # Test to set the image
        wfsImg.setImg(self.fieldXY, imageFile=self.imgFilePathIntra,
                      atype="intra")
        self.assertEqual(wfsImg.image.shape, (120, 120))
        self.assertEqual(wfsImg.fieldX, self.fieldXY[0])
        self.assertEqual(wfsImg.fieldY, self.fieldXY[1])

        # Test to update the image0
        wfsImg.updateImage0()
        self.assertEqual(np.sum(np.abs(wfsImg.image0-wfsImg.image)), 0)

        # Test to get the off axis correction
        offAxisCorrOrder = 10
        instDir = os.path.join(self.inst.instDir, self.inst.instName)
        wfsImg.getOffAxisCorr(instDir, offAxisCorrOrder)
        self.assertEqual(wfsImg.offAxisOffset, 0.001)
        self.assertEqual(wfsImg.offAxis_coeff.shape, (4, 66))
        self.assertAlmostEqual(wfsImg.offAxis_coeff[0, 0], -2.6362089*1e-3)

        # Test to make the mask list
        model = "paraxial"
        masklist = wfsImg.makeMaskList(self.inst, model)
        masklistAns = np.array([[0, 0, 1, 1], [0, 0, 0.61, 0]])
        self.assertEqual(np.sum(np.abs(masklist-masklistAns)), 0)

        model = "offAxis"
        masklist = wfsImg.makeMaskList(self.inst, model)
        masklistAns = np.array([[0, 0, 1, 1], [0, 0, 0.61, 0],
                                [-0.21240585, -0.21240585, 1.2300922, 1],
                                [-0.08784336, -0.08784336, 0.55802573, 0]])
        self.assertAlmostEqual(np.sum(np.abs(masklist-masklistAns)), 0)

        # Test the fuction to make the mask

        # The unit of boundary_thickness (boundaryT) is pixel
        boundaryT = 8
        maskScalingFactorLocal = 1
        model = "offAxis"
        wfsImg.makeMask(self.inst, model, boundaryT, maskScalingFactorLocal)
        self.assertEqual(wfsImg.pMask.shape, wfsImg.image.shape)
        self.assertEqual(wfsImg.cMask.shape, wfsImg.image.shape)
        self.assertEqual(np.sum(np.abs(wfsImg.cMask - wfsImg.pMask)), 3001)

        # Test the function of image co-center
        wfsImg.imageCoCenter(self.inst)
        xc, yc = wfsImg.getCenterAndR_ef()[0:2]
        self.assertEqual(int(xc), 63)
        self.assertEqual(int(yc), 63)

    def testFuncCompensation(self):

        # Generate a fake algorithm class
        algo = TempAlgo()
        algo.parameter = {"numTerms": 22, "offAxisPolyOrder": 10, "zobsR": 0.61}

        # Test the function of image compensation
        boundaryT = 8
        offAxisCorrOrder = 10
        zcCol = np.zeros(22)
        zcCol[3:] = self.zcCol*1e-9
        instDir = os.path.join(self.inst.instDir, self.inst.instName)

        wfsImgIntra = CompensationImageDecorator()
        wfsImgExtra = CompensationImageDecorator()
        wfsImgIntra.setImg(self.fieldXY, imageFile=self.imgFilePathIntra,
                           atype="intra")
        wfsImgExtra.setImg(self.fieldXY, imageFile=self.imgFilePathExtra,
                           atype="extra")

        for wfsImg in [wfsImgIntra, wfsImgExtra]:
            wfsImg.makeMask(self.inst, self.opticalModel, boundaryT, 1)
            wfsImg.getOffAxisCorr(instDir, offAxisCorrOrder)
            wfsImg.imageCoCenter(self.inst)
            wfsImg.compensate(self.inst, algo, zcCol, self.opticalModel)

        # Get the common region
        binaryImgIntra = wfsImgIntra.getCenterAndR_ef()[3]
        binaryImgExtra = wfsImgExtra.getCenterAndR_ef()[3]
        binaryImg = binaryImgIntra + binaryImgExtra
        binaryImg[binaryImg < 2] = 0
        binaryImg = binaryImg / 2

        # Calculate the difference
        intraImg = wfsImgIntra.getImg()
        extraImg = wfsImgExtra.getImg()
        res = np.sum(np.abs(intraImg - extraImg) * binaryImg)
        self.assertLess(res, 500)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
