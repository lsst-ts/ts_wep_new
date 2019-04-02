import os
import shutil
import unittest

from lsst.ts.wep.CamIsrWrapper import CamIsrWrapper
from lsst.ts.wep.CamDataCollector import CamDataCollector
from lsst.ts.wep.Utility import getModulePath, runProgram


class TestCamIsrWrapper(unittest.TestCase):
    """Test the CamIsrWrapper class."""

    def setUp(self):

        self.dataDir = os.path.join(getModulePath(), "tests", "tmpIsr")
        self.isrDir = os.path.join(self.dataDir, "input")
        self._makeDir(self.isrDir)

        self.camIsrWrapper = CamIsrWrapper(self.isrDir)

    def _makeDir(self, directory):

        if (not os.path.exists(directory)):
            os.makedirs(directory)

    def tearDown(self):

        shutil.rmtree(self.dataDir)

    def testConfig(self):

        fileName = self._doIsrConfig()
        isrConfigfilePath = os.path.join(self.isrDir, fileName)

        self.assertEqual(self.camIsrWrapper.doFlat, True)
        self.assertTrue(os.path.isfile(isrConfigfilePath))

        numOfLine = self._getNumOfLineInFile(isrConfigfilePath)
        self.assertEqual(numOfLine, 5)

    def _doIsrConfig(self):

        fileName = "isr_config.py"
        self.camIsrWrapper.config(doFlat=True, fileName=fileName)

        return fileName

    def _getNumOfLineInFile(self, filePath):

        with open(filePath, "r") as file:
            return sum(1 for line in file.readlines())

    def testDoIsr(self):

        # Generate the camera mapper
        camDataCollector = CamDataCollector(self.isrDir)
        camDataCollector.genPhoSimMapper()

        # Generate the fake flat images
        fakeFlatDir = os.path.join(self.dataDir, "fake_flats")
        self._makeDir(fakeFlatDir)

        detector = "R00_S22"
        self._genFakeFlat(fakeFlatDir, detector)

        # Ingest the calibration images
        calibFiles = os.path.join(fakeFlatDir, "*")
        camDataCollector.ingestCalibs(calibFiles)

        # Ingest the raw images
        imgFiles = os.path.join(getModulePath(), "tests", "testData",
                                "repackagedFiles",
                                "lsst_a_20_f5_R00_S22_E000.fits")
        camDataCollector.ingestImages(imgFiles)

        # Do the ISR configuration
        self._doIsrConfig()

        # Do the ISR
        rerunName = "run1"
        self.camIsrWrapper.doISR(self.isrDir, rerunName=rerunName)

        # Check the condition
        postIsrCcdDir = os.path.join(self.isrDir, "rerun", rerunName,
                                     "postISRCCD")
        self.assertTrue(os.path.exists(postIsrCcdDir))

    def _genFakeFlat(self, fakeFlatDir, detector):

        currWorkDir = self._getCurrWorkDir()

        self._changeWorkDir(fakeFlatDir)
        self._makeFakeFlat(detector)
        self._changeWorkDir(currWorkDir)

    def _getCurrWorkDir(self):

        return os.getcwd()

    def _changeWorkDir(self, dirPath):

        os.chdir(dirPath)

    def _makeFakeFlat(self, detector):

        command = "makeGainImages.py"
        argstring = "--detector_list %s" % detector
        runProgram(command, argstring=argstring)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
