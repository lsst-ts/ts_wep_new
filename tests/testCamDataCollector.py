import os
import shutil
import unittest

from lsst.ts.wep.CamDataCollector import CamDataCollector
from lsst.ts.wep.Utility import getModulePath, runProgram


class TestCamDataCollector(unittest.TestCase):
    """Test the CamIsrWrapper class."""

    def setUp(self):

        self.dataDir = os.path.join(getModulePath(), "tests", "tmp")
        self.isrDir = os.path.join(self.dataDir, "input")
        self._makeDir(self.isrDir)

        self.camDataCollector = CamDataCollector(self.isrDir)

    def _makeDir(self, directory):

        if (not os.path.exists(directory)):
            os.makedirs(directory)

    def tearDown(self):

        shutil.rmtree(self.dataDir)

    def testGenCamMapper(self):

        self._genMapper()

        mapperFilePath = os.path.join(self.isrDir, "_mapper")
        self.assertTrue(os.path.isfile(mapperFilePath))

        numOfLine = self._getNumOfLineInFile(mapperFilePath)
        self.assertEqual(numOfLine, 1)

    def _genMapper(self):

        self.camDataCollector.genPhoSimMapper()

    def _getNumOfLineInFile(self, filePath):

        with open(filePath, "r") as file:
            return sum(1 for line in file.readlines())

    def testIngestCalibs(self):

        # Make fake gain images
        fakeFlatDir = os.path.join(self.dataDir, "fake_flats")
        self._makeDir(fakeFlatDir)

        detector = "R00_S22"
        self._genFakeFlat(fakeFlatDir, detector)

        # Generate the mapper
        self._genMapper()

        # Do the ingestion
        calibFiles = os.path.join(fakeFlatDir, "*")
        self.camDataCollector.ingestCalibs(calibFiles)

        # Check the ingested calibration products
        calibRegistryFilePath = os.path.join(self.isrDir,
                                             "calibRegistry.sqlite3")
        self.assertTrue(os.path.exists(calibRegistryFilePath))

        flatDir = os.path.join(self.isrDir, "flat")
        self.assertTrue(os.path.exists(flatDir))

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

    def testIngestImages(self):

        self._genMapper()

        imgFiles = os.path.join(getModulePath(), "tests", "testData",
                                "repackagedFiles",
                                "lsst_a_20_f5_R00_S22_E000.fits")
        self.camDataCollector.ingestImages(imgFiles)

        # Check the ingested calibration products
        registryFilePath = os.path.join(self.isrDir, "registry.sqlite3")
        self.assertTrue(os.path.exists(registryFilePath))

        rawDir = os.path.join(self.isrDir, "raw")
        self.assertTrue(os.path.exists(rawDir))


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
