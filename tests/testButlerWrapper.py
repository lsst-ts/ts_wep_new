import os
import unittest

from lsst.ts.wep.ButlerWrapper import ButlerWrapper
from lsst.ts.wep.Utility import getModulePath


class TestButlerWrapper(unittest.TestCase):
    """Test the butler wrapper class."""

    def setUp(self):

        self.inputs = os.path.join(getModulePath(), "tests", "testData",
                                   "repackagedPhoSimData")
        self.butlerWrapper = ButlerWrapper(inputs=self.inputs)

    def testGetPhoSimRpkgRawExp(self):

        visit, raft, sensor = self._getDefaultSurveyMetaData()
        exposure = self.butlerWrapper.getPhoSimRpkgRawExp(visit, raft, sensor)
        self.assertEqual(exposure.getDimensions()[0], 4176)
        self.assertEqual(exposure.getDimensions()[1], 4020)

    def _getDefaultSurveyMetaData(self):

        visit = 20
        raft = "R00"
        sensor = "S22"

        return visit, raft, sensor

    def testGetPhoSimPostIsrCcd(self):

        postIsrCcd = self._getPhoSimPostIsrCcd()
        self.assertEqual(postIsrCcd.getDimensions()[0], 4072)
        self.assertEqual(postIsrCcd.getDimensions()[1], 4000)

    def _getPhoSimPostIsrCcd(self):

        visit, raft, sensor = self._getDefaultSurveyMetaData()
        postIsrCcd = self.butlerWrapper.getPhoSimPostIsrCcd(visit, raft, sensor)

        return postIsrCcd

    def testSetInputsAndOutputs(self):

        self.butlerWrapper.setInputsAndOutputs(inputs=self.inputs)

        postIsrCcd = self._getPhoSimPostIsrCcd()
        self.assertEqual(postIsrCcd.getDimensions()[0], 4072)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
