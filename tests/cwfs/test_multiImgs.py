# This file is part of ts_wep.
#
# Developed for the LSST Telescope and Site Systems.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os
import time
import numpy as np
import unittest

from lsst.ts.wep.cwfs.Instrument import Instrument
from lsst.ts.wep.cwfs.Algorithm import Algorithm
from lsst.ts.wep.cwfs.CompensableImage import CompensableImage
from lsst.ts.wep.PlotUtil import plotImage
from lsst.ts.wep.Utility import getModulePath, getConfigDir, DefocalType, CamType


class WepFile(object):
    def __init__(
        self,
        imageFolderName,
        imageName,
        fieldXY,
        useAlgorithm,
        orientation,
        validationPath,
    ):

        self.imageFolderName = imageFolderName
        self.imageName = imageName
        self.fieldXY = fieldXY
        self.useAlgorithm = useAlgorithm
        self.orientation = orientation

        # Get the refFilePath
        refFileName = imageFolderName + "_" + imageName + useAlgorithm + ".txt"
        self.refFilePath = os.path.join(validationPath, refFileName)


class DataWep(object):
    def __init__(self, instFolder, algoFolderPath, imageFolder, wepFile):

        self.instFolder = instFolder
        self.algoFolderPath = algoFolderPath
        self.useAlgorithm = wepFile.useAlgorithm
        self.imageFolderPath = os.path.join(imageFolder, wepFile.imageFolderName)
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
        self.imageFolderPath = os.path.join(
            modulePath, "tests", "testData", "testImages"
        )

        # Set the tolerance
        self.tor = 4

        # Restart time
        self.startTime = time.time()
        self.difference = 0
        self.validationDir = os.path.join(
            modulePath, "tests", "testData", "testImages", "validation", "simulation"
        )

    def tearDown(self):

        # Calculate the time of test case
        t = time.time() - self.startTime
        print("%s: %.3f s. Differece is %.3f." % (self.id(), t, self.difference))

    def _generateTestCase(
        self,
        imageFolderName,
        imageName,
        fieldXY,
        useAlgorithm,
        orientation,
        validationPath,
    ):

        caseTest = WepFile(
            imageFolderName,
            imageName,
            fieldXY,
            useAlgorithm,
            orientation,
            validationPath,
        )
        case = DataWep(
            self.instFolder, self.algoFolderPath, self.imageFolderPath, caseTest
        )

        return case

    def _compareCalculation(self, dataWEP, tor):

        # Run WEP to get Zk
        zer4UpNm = self._runWEP(
            dataWEP.instFolder,
            dataWEP.algoFolderPath,
            dataWEP.useAlgorithm,
            dataWEP.imageFolderPath,
            dataWEP.intra_image_name,
            dataWEP.extra_image_name,
            dataWEP.fieldXY,
            dataWEP.orientation,
        )

        # Load the reference data
        refZer4UpNm = np.loadtxt(dataWEP.refFilePath)

        # Compare the result
        difference = np.max(np.abs(zer4UpNm - refZer4UpNm))

        if difference <= tor:
            result = "true"
        else:
            result = "false"

        self.difference = difference

        return result

    def _runWEP(
        self,
        instDir,
        algoFolderPath,
        useAlgorithm,
        imageFolderPath,
        intra_image_name,
        extra_image_name,
        fieldXY,
        opticalModel,
        showFig=False,
    ):

        # Image files Path
        intra_image_file = os.path.join(imageFolderPath, intra_image_name)
        extra_image_file = os.path.join(imageFolderPath, extra_image_name)

        # There is the difference between intra and extra images
        # I1: intra_focal images, I2: extra_focal Images
        I1 = CompensableImage()
        I2 = CompensableImage()

        I1.setImg(fieldXY, DefocalType.Intra, imageFile=intra_image_file)
        I2.setImg(fieldXY, DefocalType.Extra, imageFile=extra_image_file)

        # Set the instrument
        inst = Instrument(instDir)
        inst.config(CamType.LsstCam, I1.getImgSizeInPix(), announcedDefocalDisInMm=1.0)

        # Define the algorithm to be used.
        algo = Algorithm(algoFolderPath)
        algo.config(useAlgorithm, inst, debugLevel=0)

        # Plot the original wavefront images
        if showFig:
            plotImage(I1.image, title="intra image")
            plotImage(I2.image, title="extra image")

        # Run it
        algo.runIt(I1, I2, opticalModel, tol=1e-3)

        # Show the Zernikes Zn (n>=4)
        algo.outZer4Up(showPlot=False)

        # Plot the final conservated images and wavefront
        if showFig:
            plotImage(I1.image, title="Compensated intra image")
            plotImage(I2.image, title="Compensated extra image")

            # Plot the Wavefront
            plotImage(algo.wcomp, title="Final wavefront")
            plotImage(
                algo.wcomp,
                title="Final wavefront with pupil mask applied",
                mask=algo.pMask,
            )

        # Return the Zernikes Zn (n>=4)
        return algo.getZer4UpInNm()

    def testCase1(self):

        case = self._generateTestCase(
            "LSST_NE_SN25",
            "z11_0.25_",
            (1.185, 1.185),
            "exp",
            "offAxis",
            self.validationDir,
        )
        result = self._compareCalculation(case, self.tor)
        self.assertEqual(result, "true")

    def testCase2(self):

        case = self._generateTestCase(
            "LSST_NE_SN25",
            "z11_0.25_",
            (1.185, 1.185),
            "fft",
            "offAxis",
            self.validationDir,
        )
        result = self._compareCalculation(case, self.tor)
        self.assertEqual(result, "true")

    def testCase3(self):

        case = self._generateTestCase(
            "F1.23_1mm_v61", "z7_0.25_", (0, 0), "fft", "paraxial", self.validationDir
        )
        result = self._compareCalculation(case, self.tor)
        self.assertEqual(result, "true")

    def testCase4(self):

        case = self._generateTestCase(
            "LSST_C_SN26", "z7_0.25_", (0, 0), "fft", "onAxis", self.validationDir
        )
        result = self._compareCalculation(case, self.tor)
        self.assertEqual(result, "true")

    def testCase5(self):

        case = self._generateTestCase(
            "LSST_C_SN26", "z7_0.25_", (0, 0), "exp", "onAxis", self.validationDir
        )
        result = self._compareCalculation(case, self.tor)
        self.assertEqual(result, "true")


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
