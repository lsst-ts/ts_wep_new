import unittest

from lsst.ts.wep.bsc.LsstCam import LsstCam


class TestLsstCam(unittest.TestCase):
    """Test the LsstCam class."""

    def setUp(self):

        self.cam = LsstCam()
        self.cam.setObsMetaData(0, 0, 0, mjd=59580.0)

    def testGetWfsCcdList(self):

        wfsList = self.cam.getWfsCcdList()
        self.assertEqual(len(wfsList), 8)

    def testGetWavefrontSensor(self):

        wfsData = self.cam.getWavefrontSensor()
        self.assertEqual(len(wfsData), 8)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
