import unittest

from lsst.ts.wep.bsc.LsstFamCam import LsstFamCam


class TestLsstFamCam(unittest.TestCase):
    """Test the LsstCam class."""

    def setUp(self):

        self.cam = LsstFamCam()
        self.cam.setObsMetaData(0, 0, 0, mjd=59580.0)

    def testGetWfsCcdList(self):

        wfsList = self.cam.getWfsCcdList()
        self.assertEqual(len(wfsList), 189)

    def testGetWavefrontSensor(self):

        wfsData = self.cam.getWavefrontSensor()
        self.assertEqual(len(wfsData), 189)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
