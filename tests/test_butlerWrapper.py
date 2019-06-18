import os
import numpy as np
import unittest

from lsst.ts.wep.ButlerWrapper import ButlerWrapper
from lsst.ts.wep.Utility import getModulePath


class TestButlerWrapper(unittest.TestCase):
    """Test the butler wrapper class."""

    def setUp(self):

        self.inputs = os.path.join(getModulePath(), "tests", "testData",
                                   "repackagedPhoSimData")
        self.butlerWrapper = ButlerWrapper(self.inputs)

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

        self.butlerWrapper.setInputsAndOutputs(inputs=self.inputs)

        exposure = self._getRawExp()
        self.assertEqual(exposure.getDimensions().getX(), 4176)

    def testGetImageData(self):

        exposure = self._getRawExp()

        image = ButlerWrapper.getImageData(exposure)
        self.assertTrue(isinstance(image, np.ndarray))
        self.assertEqual(image.shape, (4020, 4176))

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
