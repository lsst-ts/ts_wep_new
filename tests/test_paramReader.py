import os
import shutil
import numpy as np
import unittest

from lsst.ts.wep.ParamReader import ParamReader
from lsst.ts.wep.Utility import getModulePath


class TestParamReader(unittest.TestCase):
    """Test the ParamReaderYaml class."""

    def setUp(self):

        testDir = os.path.join(getModulePath(), "tests")
        self.configDir = os.path.join(testDir, "testData")
        self.fileName = "testConfigFile.yaml"

        filePath = os.path.join(self.configDir, self.fileName)
        self.paramReader = ParamReader(filePath=filePath)

        self.testTempDir = os.path.join(testDir, "tmp")
        self._makeDir(self.testTempDir)

    def _makeDir(self, directory):

        if (not os.path.exists(directory)):
            os.makedirs(directory)

    def tearDown(self):

        shutil.rmtree(self.testTempDir)

    def testGetSetting(self):

        znmax = self.paramReader.getSetting("znmax")
        self.assertEqual(znmax, 22)

    def testGetFilePath(self):

        ansFilePath = os.path.join(self.configDir, self.fileName)
        self.assertEqual(self.paramReader.getFilePath(), ansFilePath)

    def testSetFilePath(self):

        fileName = "test.yaml"
        filePath = os.path.join(self.configDir, fileName)
        self.paramReader.setFilePath(filePath)

        self.assertEqual(self.paramReader.getFilePath(), filePath)

    def testGetContent(self):

        content = self.paramReader.getContent()
        self.assertTrue(isinstance(content, dict))

    def testGetContentWithDefaultSetting(self):

        paramReader = ParamReader()

        content = paramReader.getContent()
        self.assertTrue(isinstance(content, dict))

    def testWriteMatToFile(self):

        self._writeMatToFile()

        numOfFile = self._getNumOfFileInFolder(self.testTempDir)
        self.assertEqual(numOfFile, 1)

    def _writeMatToFile(self):

        mat = np.random.rand(3, 4, 5)
        filePath = os.path.join(self.testTempDir, "temp.yaml")
        ParamReader.writeMatToFile(mat, filePath)

        return mat, filePath

    def _getNumOfFileInFolder(self, folder):

        return len([name for name in os.listdir(folder)
                   if os.path.isfile(os.path.join(folder, name))])

    def testWriteMatToFileWithWrongFileFormat(self):

        wrongFilePath = os.path.join(self.testTempDir, "temp.txt")
        self.assertRaises(ValueError, ParamReader.writeMatToFile,
                          np.ones(4), wrongFilePath)

    def testGetMatContent(self):

        mat, filePath = self._writeMatToFile()

        self.paramReader.setFilePath(filePath)
        matInYamlFile = self.paramReader.getMatContent()

        delta = np.sum(np.abs(matInYamlFile - mat))
        self.assertLess(delta, 1e-10)

    def testGetMatContentWithDefaultSetting(self):

        paramReader = ParamReader()
        matInYamlFile = paramReader.getMatContent()

        self.assertTrue(isinstance(matInYamlFile, np.ndarray))
        self.assertEqual(len(matInYamlFile), 0)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
