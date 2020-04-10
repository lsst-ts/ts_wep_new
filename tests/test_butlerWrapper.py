import os
import numpy as np
import tempfile
import unittest

from lsst.ts.wep.CamDataCollector import CamDataCollector
from lsst.ts.wep.ButlerWrapper import ButlerWrapper
from lsst.ts.wep.Utility import getModulePath


class TestButlerWrapper(unittest.TestCase):
    """Test the butler wrapper class."""

    @classmethod
    def setUpClass(cls):

        dataDirPath = os.path.join(getModulePath(), "tests")
        cls.dataDir = tempfile.TemporaryDirectory(dir=dataDirPath)
        cls._ingestImages()

        cls.butlerWrapper = ButlerWrapper(cls.dataDir.name)

    @classmethod
    def _ingestImages(cls):

        # Generate the camera mapper
        camDataCollector = CamDataCollector(cls.dataDir.name)
        camDataCollector.genPhoSimMapper()

        # Ingest the E image
        imgFilesEimg = os.path.join(getModulePath(), "tests", "testData",
                                    "repackagedFiles",
                                    "lsst_e_9006001_f1_R22_S00_E000.fits.gz")
        camDataCollector.ingestEimages(imgFilesEimg)

        # Ingest the amplifier image
        imgFilesRaw = os.path.join(getModulePath(), "tests", "testData",
                                   "repackagedFiles",
                                   "lsst_a_20_f5_R00_S22_E000.fits")
        camDataCollector.ingestImages(imgFilesRaw)

    @classmethod
    def tearDownClass(cls):

        cls.dataDir.cleanup()

    def testGetRawExp(self):

        exposure = self._getRawExp()
        self.assertEqual(exposure.getDimensions().getX(), 4176)
        self.assertEqual(exposure.getDimensions().getY(), 4020)

    def _getRawExp(self):

        visit = 20
        raft = "R00"
        sensor = "S22"
        return self.butlerWrapper.getRawExp(visit, raft, sensor)

    def testGetEimage(self):

        exposure = self._getEimage()
        self.assertEqual(exposure.getDimensions().getX(), 4072)
        self.assertEqual(exposure.getDimensions().getY(), 4000)

    def _getEimage(self):

        visit = 9006001
        raft = "R22"
        sensor = "S00"
        return self.butlerWrapper.getEimage(visit, raft, sensor)

    def testSetInputsAndOutputs(self):

        self.butlerWrapper.setInputsAndOutputs(self.dataDir.name)

        eimg = self._getEimage()
        self.assertEqual(eimg.getDimensions().getX(), 4072)

    def testGetImageData(self):

        eimg = self._getEimage()

        image = ButlerWrapper.getImageData(eimg)
        self.assertTrue(isinstance(image, np.ndarray))
        self.assertEqual(image.shape, (4000, 4072))

    def testExtendDataId(self):

        dataId = dict()

        snap = 0
        aFilter = "u"
        self.butlerWrapper._extendDataId(dataId, snap=snap, aFilter=aFilter)

        self.assertEqual(dataId["snap"], snap)
        self.assertEqual(dataId["filter"], aFilter)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
