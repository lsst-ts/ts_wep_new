import os
import sys
import time
import numpy as np
import unittest

from lsst.ts.wep.cwfs import Instrument, Algorithm, CompensationImageDecorator
from lsst.ts.wep.cwfs.Tool import plotImage
from lsst.ts.wep.Utility import getModulePath


def runWEP(instruFolder, algoFolderPath, instruName, useAlgorithm,
           imageFolderPath, intra_image_name, extra_image_name,
           fieldXY, opticalModel, showFig=False, showConf=False,
           filename=None):
    """Calculate the coefficients of normal/ annular Zernike polynomials based
    on the provided instrument, algorithm, and optical model.

    Parameters
    ----------
    instruFolder : str
        Path to instrument folder.
    algoFolderPath : str
        Path to algorithm folder.
    instruName : str
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
    showConf : bool, optional
        Decide to show the configuration or not. (the default is False.)
    filename : str, optional
        Name of output file. (the default is None.)

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
    I1 = CompensationImageDecorator.CompensationImageDecorator()
    I2 = CompensationImageDecorator.CompensationImageDecorator()

    I1.setImg(fieldXY, imageFile=intra_image_file, atype="intra")
    I2.setImg(fieldXY, imageFile=extra_image_file, atype="extra")

    # Set the instrument
    inst = Instrument.Instrument(instruFolder)
    inst.config(instruName, I1.sizeinPix)

    # Define the algorithm to be used.
    algo = Algorithm.Algorithm(algoFolderPath)
    algo.config(useAlgorithm, inst, debugLevel=0)

    # Plot the original wavefront images
    if (showFig):
        plotImage(I1.image, title="intra image")
        plotImage(I2.image, title="extra image")

    # Run it
    algo.runIt(inst, I1, I2, opticalModel, tol=1e-3)

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

    # Show the information of image, algorithm, and instrument
    if (showConf):
        _outParam(algo, inst, I1, I2, opticalModel, filename)

    # Return the Zernikes Zn (n>=4)
    return algo.zer4UpNm


def _outParam(algo, inst, I1, I2, opticalModel, filename=None):
    """Put the information of images, instrument, and algorithm on terminal or
    file.

    Parameters
    ----------
    algo : Algorithm
        Algorithm to solve the Poisson's equation in the transport of intensity
        equation (TIE).
    inst : Instrument
        Instrument to use.
    I1 : Image
        Intra- or extra-focal image.
    I2 : Image
        Intra- or extra-focal image.
    opticalModel : str
        Optical model. It can be "paraxial", "onAxis", or "offAxis".
    filename : str, optional
        Name of output file. (the default is None.)
    """

    # Write the parameters into a file if needed.
    if (filename is not None):
        fout = open(filename, "w")
    else:
        fout = sys.stdout

    # Write the information of image and optical model
    fout.write("intra image: \t %s \t field in deg =(%6.3f, %6.3f)\n" %
               (I1.name, I1.fieldX, I1.fieldY))
    fout.write("extra image: \t %s \t field in deg =(%6.3f, %6.3f)\n" %
               (I2.name, I2.fieldX, I2.fieldY))
    fout.write("Using optical model:\t %s\n" % opticalModel)

    # Read the instrument file
    _readConfigFile(fout, inst, "instrument")

    # Read the algorithm file
    _readConfigFile(fout, algo, "algorithm")

    # Close the file
    if (filename is not None):
        fout.close()


def _readConfigFile(fout, config, configName):
    """Read the configuration file

    Parameters
    ----------
    fout : file
        File instance.
    config : Instrument or Algorithm
        Instance of configuration.
    configName : str
        Name of configuration.
    """

    # Create a new line
    fout.write("\n")

    # Open the file
    fconfig = open(config.filename)
    fout.write("---" + configName + " file: --- %s ----------\n" % config.filename)

    # Read the file information
    iscomment = False
    for line in fconfig:
        line = line.strip()
        if (line.startswith("###")):
            iscomment = ~iscomment
        if (not(line.startswith("#")) and (not iscomment) and len(line) > 0):
            fout.write(line + "\n")

    # Close the file
    fconfig.close()


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
        self.instFolder = os.path.join(modulePath, "configData", "cwfs",
                                       "instruData")
        self.algoFolderPath = os.path.join(modulePath, "configData", "cwfs",
                                           "algo")
        self.imageFolderPath = os.path.join(modulePath, "tests", "testData",
                                            "testImages")
        self.instName = "lsst"

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
