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
        self.assertEqual(exposure.getDimensions()[0], 4176)
        self.assertEqual(exposure.getDimensions()[1], 4020)

    def _getRawExp(self):

        visit, raft, sensor = self._getDefaultSurveyMetaData()
        exposure = self.butlerWrapper.getRawExp(visit, raft, sensor)

        return exposure

    def _getDefaultSurveyMetaData(self):

        visit = 20
        raft = "R00"
        sensor = "S22"

        return visit, raft, sensor

    def testSetInputsAndOutputs(self):

        self.butlerWrapper.setInputsAndOutputs(inputs=self.inputs)

        exposure = self._getRawExp()
        self.assertEqual(exposure.getDimensions()[0], 4176)

    def testGetImageData(self):

        exposure = self._getRawExp()

        image = ButlerWrapper.getImageData(exposure)
        self.assertTrue(isinstance(image, np.ndarray))
        self.assertEqual(image.shape, (4020, 4176))


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
