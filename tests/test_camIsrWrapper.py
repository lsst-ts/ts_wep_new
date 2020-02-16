import os
import shutil
import unittest

from lsst.ts.wep.CamIsrWrapper import CamIsrWrapper
from lsst.ts.wep.CamDataCollector import CamDataCollector
from lsst.ts.wep.Utility import getModulePath, runProgram


class TestCamIsrWrapper(unittest.TestCase):
    """Test the CamIsrWrapper class."""

    def setUp(self):

        testDir = os.path.join(getModulePath(), "tests")
        self.dataDir = os.path.join(testDir, "tmpIsr")
        self.isrDir = os.path.join(self.dataDir, "input")
        self._makeDir(self.isrDir)

        self.repackagedTestData = os.path.join(testDir, "testData",
                                               "repackagedFiles")

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

        # Get the camDataCollector and ingest the calibs
        detector = "R00_S22"
        camDataCollector = self._getCamDataCollectorAndIngestCalibs(
            detector)

        # Ingest the raw images and do the ISR
        imgFileName = "lsst_a_20_f5_R00_S22_E000.fits"
        rerunName = "run1"
        self._ingestRawAndDoIsr(imgFileName, camDataCollector, rerunName,
                                doIsrConfig=True)

        # Check the condition
        postIsrCcdDir = os.path.join(self.isrDir, "rerun", rerunName,
                                     "postISRCCD")
        self.assertTrue(os.path.exists(postIsrCcdDir))

        numOfDir = self._getNumOfDir(postIsrCcdDir)
        self.assertEqual(numOfDir, 1)

    def _getCamDataCollectorAndIngestCalibs(self, detector):

        # Generate the camera mapper
        camDataCollector = CamDataCollector(self.isrDir)
        camDataCollector.genPhoSimMapper()

        # Generate the fake flat images
        fakeFlatDir = os.path.join(self.dataDir, "fake_flats")
        self._makeDir(fakeFlatDir)

        self._genFakeFlat(fakeFlatDir, detector)

        # Ingest the calibration images
        calibFiles = os.path.join(fakeFlatDir, "*")
        camDataCollector.ingestCalibs(calibFiles)

        return camDataCollector

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

    def _ingestRawAndDoIsr(self, imgFileName, camDataCollector, rerunName,
                           doIsrConfig=False):

        imgFiles = os.path.join(self.repackagedTestData, imgFileName)
        camDataCollector.ingestImages(imgFiles)

        # Do the ISR configuration
        if (doIsrConfig is True):
            self._doIsrConfig()

        # Do the ISR
        self.camIsrWrapper.doISR(self.isrDir, rerunName=rerunName)

    def _getNumOfDir(self, dirPath):

        numOfDir = sum(os.path.isdir(os.path.join(dirPath, aDir))
                       for aDir in os.listdir(path=dirPath))

        return numOfDir

    @unittest.skip
    def testDoIsrContinuous(self):

        # In lsst_distrib w_2020_06
        # There is the bug in upstream that the user can not do the ISR
        # continuous. Report this bug to DM team and wait for the reply.

        # Get the camDataCollector and ingest the calibs
        detector = "R00_S22 R22_S10"
        camDataCollector = self._getCamDataCollectorAndIngestCalibs(
            detector)

        # Ingest the first raw images and do the ISR
        imgFileName = "lsst_a_20_f5_R00_S22_E000.fits"
        rerunName = "run1"
        self._ingestRawAndDoIsr(imgFileName, camDataCollector, rerunName,
                                doIsrConfig=True)

        # Check the condition
        postIsrCcdDir = os.path.join(self.isrDir, "rerun", rerunName,
                                     "postISRCCD")
        numOfDir = self._getNumOfDir(postIsrCcdDir)
        self.assertEqual(numOfDir, 1)

        # Ingest the second raw image and do the ISR again
        imgFileName = "lsst_a_9005000_f1_R22_S10_E000.fits"
        self._ingestRawAndDoIsr(imgFileName, camDataCollector, rerunName,
                                doIsrConfig=False)

        # Check the condition again
        numOfDir = self._getNumOfDir(postIsrCcdDir)
        self.assertEqual(numOfDir, 2)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
