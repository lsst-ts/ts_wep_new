import unittest

from lsst.ts.wep.Utility import CamType
from lsst.ts.wep.bsc.CamFactory import CamFactory
from lsst.ts.wep.bsc.LsstCam import LsstCam
from lsst.ts.wep.bsc.LsstFamCam import LsstFamCam
from lsst.ts.wep.bsc.ComCam import ComCam


class TestCamFactory(unittest.TestCase):
    """Test the CamFactory class."""

    def testCreateCam(self):

        lsstCam = CamFactory.createCam(CamType.LsstCam)
        self.assertTrue(isinstance(lsstCam, LsstCam))

        lsstFamCam = CamFactory.createCam(CamType.LsstFamCam)
        self.assertTrue(isinstance(lsstFamCam, LsstFamCam))

        comCam = CamFactory.createCam(CamType.ComCam)
        self.assertTrue(isinstance(comCam, ComCam))

        self.assertRaises(ValueError, CamFactory.createCam, "wrongType")


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
