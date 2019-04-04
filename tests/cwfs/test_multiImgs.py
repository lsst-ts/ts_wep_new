import os
import time
import numpy as np
import unittest

from lsst.ts.wep.cwfs.Instrument import Instrument
from lsst.ts.wep.cwfs.Algorithm import Algorithm
from lsst.ts.wep.cwfs.CompensationImageDecorator import CompensationImageDecorator
from lsst.ts.wep.cwfs.Tool import plotImage
from lsst.ts.wep.Utility import getModulePath, getConfigDir, DefocalType


def runWEP(instDir, algoFolderPath, instName, useAlgorithm,
           imageFolderPath, intra_image_name, extra_image_name,
           fieldXY, opticalModel, showFig=False):
    """Calculate the coefficients of normal/ annular Zernike polynomials based
    on the provided instrument, algorithm, and optical model.

    Parameters
    ----------
    instDir : str
        Path to instrument folder.
    algoFolderPath : str
        Path to algorithm folder.
    instName : str
        Instrument name. It is "lsst" in the baseline.
    useAlgorithm : str
        Algorithm to solve the Poisson's equation in the transport of intensity
        equation (TIE). It can be "fft" or "exp" here.
    imageFolderPath : str
        Path to image folder.
    intra_image_name : str
        File name of intra-focal image.
    extra_image_name : str
        File name of extra-focal image.
    fieldXY : tuple
        Position of donut on the focal plane in degree for intra- and
        extra-focal images.
    opticalModel : str
        Optical model. It can be "paraxial", "onAxis", or "offAxis".
    showFig : bool, optional
        Show the wavefront image and compenstated image or not. (the default is
        False.)

    Returns
    -------
    numpy.ndarray
        Coefficients of Zernike polynomials (z4 - z22).
    """

    # Image files Path
    intra_image_file = os.path.join(imageFolderPath, intra_image_name)
    extra_image_file = os.path.join(imageFolderPath, extra_image_name)

    # There is the difference between intra and extra images
    # I1: intra_focal images, I2: extra_focal Images
    I1 = CompensationImageDecorator()
    I2 = CompensationImageDecorator()

    I1.setImg(fieldXY, DefocalType.Intra, imageFile=intra_image_file)
    I2.setImg(fieldXY, DefocalType.Extra, imageFile=extra_image_file)

    # Set the instrument
    inst = Instrument(instDir)
    inst.config(instName, I1.getImgSizeInPix())

    # Define the algorithm to be used.
    algo = Algorithm(algoFolderPath)
    algo.config(useAlgorithm, inst, debugLevel=0)

    # Plot the original wavefront images
    if (showFig):
        plotImage(I1.image, title="intra image")
        plotImage(I2.image, title="extra image")

    # Run it
    algo.runIt(I1, I2, opticalModel, tol=1e-3)

    # Show the Zernikes Zn (n>=4)
    algo.outZer4Up(showPlot=False)

    # Plot the final conservated images and wavefront
    if (showFig):
        plotImage(I1.image, title="Compensated intra image")
        plotImage(I2.image, title="Compensated extra image")

        # Plot the Wavefront
        plotImage(algo.wcomp, title="Final wavefront")
        plotImage(algo.wcomp, title="Final wavefront with pupil mask applied",
                  mask=algo.pMask)

    # Return the Zernikes Zn (n>=4)
    return algo.getZer4UpInNm()


class WepFile(object):

    def __init__(self, imageFolderName, imageName, fieldXY, useAlgorithm,
                 orientation, validationPath):

        self.imageFolderName = imageFolderName
        self.imageName = imageName
        self.fieldXY = fieldXY
        self.useAlgorithm = useAlgorithm
        self.orientation = orientation

        # Get the refFilePath
        refFileName = imageFolderName + "_" + imageName + useAlgorithm + ".txt"
        self.refFilePath = os.path.join(validationPath, refFileName)


class DataWep(object):

    def __init__(self, instFolder, algoFolderPath, instName, imageFolder,
                 wepFile):

        self.instruFolder = instFolder
        self.algoFolderPath = algoFolderPath
        self.instruName = instName
        self.useAlgorithm = wepFile.useAlgorithm
        self.imageFolderPath = os.path.join(imageFolder,
                                            wepFile.imageFolderName)
        self.intra_image_name = wepFile.imageName + "intra.txt"
        self.extra_image_name = wepFile.imageName + "extra.txt"
        self.fieldXY = wepFile.fieldXY
        self.orientation = wepFile.orientation
        self.refFilePath = wepFile.refFilePath


class TestWepWithMultiImgs(unittest.TestCase):

    def setUp(self):

        # Set the path of module and the setting directories
        modulePath = getModulePath()

        cwfsConfigDir = os.path.join(getConfigDir(), "cwfs")
        self.instFolder = os.path.join(cwfsConfigDir, "instData")
        self.algoFolderPath = os.path.join(cwfsConfigDir, "algo")
        self.imageFolderPath = os.path.join(modulePath, "tests", "testData",
                                            "testImages")
        self.instName = "lsst10"

        # Set the tolerance
        self.tor = 3

        # Restart time
        self.startTime = time.time()
        self.difference = 0
        self.validationDir = os.path.join(modulePath, "tests", "testData",
                                          "testImages", "validation")

    def tearDown(self):

        # Calculate the time of test case
        t = time.time() - self.startTime
        print("%s: %.3f s. Differece is %.3f." % (self.id(), t,
              self.difference))

    def _generateTestCase(self, imageFolderName, imageName, fieldXY,
                          useAlgorithm, orientation, validationPath):

        caseTest = WepFile(imageFolderName, imageName, fieldXY, useAlgorithm,
                           orientation, validationPath)
        case = DataWep(self.instFolder, self.algoFolderPath, self.instName,
                       self.imageFolderPath, caseTest)

        return case

    def _compareCalculation(self, dataWEP, tor):

        # Run WEP to get Zk
        zer4UpNm = runWEP(dataWEP.instruFolder, dataWEP.algoFolderPath,
                          dataWEP.instruName, dataWEP.useAlgorithm,
                          dataWEP.imageFolderPath, dataWEP.intra_image_name,
                          dataWEP.extra_image_name, dataWEP.fieldXY,
                          dataWEP.orientation)

        # Load the reference data
        refZer4UpNm = np.loadtxt(dataWEP.refFilePath)

        # Compare the result
        difference = np.sum((zer4UpNm-refZer4UpNm)**2)

        if difference <= tor:
            result = "true"
        else:
            result = "false"

        self.difference = difference

        return result

    def testCase1(self):

        case = self._generateTestCase("LSST_NE_SN25", "z11_0.25_",
                                      (1.185, 1.185), "exp", "offAxis",
                                      self.validationDir)
        result = self._compareCalculation(case, self.tor)
        self.assertEqual(result, "true")

    def testCase2(self):

        case = self._generateTestCase("LSST_NE_SN25", "z11_0.25_",
                                      (1.185, 1.185), "fft", "offAxis",
                                      self.validationDir)
        result = self._compareCalculation(case, self.tor)
        self.assertEqual(result, "true")

    def testCase3(self):

        case = self._generateTestCase("F1.23_1mm_v61", "z7_0.25_", (0, 0),
                                      "fft", "paraxial", self.validationDir)
        result = self._compareCalculation(case, self.tor)
        self.assertEqual(result, "true")

    def testCase4(self):

        case = self._generateTestCase("LSST_C_SN26", "z7_0.25_", (0, 0),
                                      "fft", "onAxis", self.validationDir)
        result = self._compareCalculation(case, self.tor)
        self.assertEqual(result, "true")

    def testCase5(self):

        case = self._generateTestCase("LSST_C_SN26", "z7_0.25_", (0, 0),
                                      "exp", "onAxis", self.validationDir)
        result = self._compareCalculation(case, self.tor)
        self.assertEqual(result, "true")


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
