import os
import tempfile
import unittest

from lsst.ts.wep.CamDataCollector import CamDataCollector
from lsst.ts.wep.Utility import getModulePath, runProgram


class TestCamDataCollector(unittest.TestCase):
    """Test the CamIsrWrapper class."""

    def setUp(self):

        testDir = os.path.join(getModulePath(), "tests")
        self.dataDir = tempfile.TemporaryDirectory(dir=testDir)

        self.isrDir = tempfile.TemporaryDirectory(dir=self.dataDir.name)

        self.camDataCollector = CamDataCollector(self.isrDir.name)

    def tearDown(self):

        self.dataDir.cleanup()

    def testGenCamMapper(self):

        self._genMapper()

        mapperFilePath = os.path.join(self.isrDir.name, "_mapper")
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
        fakeFlatDir = tempfile.TemporaryDirectory(dir=self.dataDir.name)

        detector = "R00_S22"
        self._genFakeFlat(fakeFlatDir.name, detector)

        # Generate the mapper
        self._genMapper()

        # Do the ingestion
        calibFiles = os.path.join(fakeFlatDir.name, "*")
        self.camDataCollector.ingestCalibs(calibFiles)

        # Check the ingested calibration products
        calibRegistryFilePath = os.path.join(self.isrDir.name,
                                             "calibRegistry.sqlite3")
        self.assertTrue(os.path.exists(calibRegistryFilePath))

        flatDir = os.path.join(self.isrDir.name, "flat")
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

        # Check the ingested files
        self._checkIngestion("raw")

    def _checkIngestion(self, imgType):

        registryFilePath = os.path.join(self.isrDir.name, "registry.sqlite3")
        self.assertTrue(os.path.exists(registryFilePath))

        rawDir = os.path.join(self.isrDir.name, imgType)
        self.assertTrue(os.path.exists(rawDir))

    def testIngestEimages(self):

        self._genMapper()

        imgFiles = os.path.join(getModulePath(), "tests", "testData",
                                "repackagedFiles",
                                "lsst_e_9006001_f1_R22_S00_E000.fits.gz")
        self.camDataCollector.ingestEimages(imgFiles)

        # Check the ingested files
        self._checkIngestion("eimage")


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
