import os
import shutil
import unittest

from lsst.ts.wep.CamIsrWrapper import CamIsrWrapper
from lsst.ts.wep.Utility import getModulePath


class TestCamIsrWrapper(unittest.TestCase):
    """Test the CamIsrWrapper class."""

    def setUp(self):
        
        self.dataDir = os.path.join(getModulePath(), "tests", "tmp")
        self.isrDir = os.path.join(self.dataDir, "input")
        self._makeDir(self.isrDir)

        self.camIsrWrapper = CamIsrWrapper(destDir=self.isrDir)

    def _makeDir(self, directory):

        if (not os.path.exists(directory)):
            os.makedirs(directory)

    def tearDown(self):

        shutil.rmtree(self.dataDir)

    def testConfig(self):

        fileName="isr_config.py"
        self.camIsrWrapper.config(doFlat=True, fileName=fileName)

        isrConfigfilePath = os.path.join(self.isrDir, fileName)

        self.assertEqual(self.camIsrWrapper.doFlat, True)
        self.assertTrue(os.path.isfile(isrConfigfilePath))

        numOfLine = self._getNumOfLineInFile(isrConfigfilePath)
        self.assertEqual(numOfLine, 5)

    def _getNumOfLineInFile(self, filePath):

        with open(filePath, "r") as file:
            return sum(1 for line in file.readlines())

    def testDoIsr(self):
        pass


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
