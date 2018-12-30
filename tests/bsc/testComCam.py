import unittest

from lsst.ts.wep.bsc.StarData import StarData
from lsst.ts.wep.bsc.ComCam import ComCam
from lsst.ts.wep.Utility import FilterType



class TestComCam(unittest.TestCase):
    """Test the ComCam class."""

    camera = None

    # Stars
    stars = None

    def setUp(self):

        # Boresight (unit: degree)
        ra = 0.0    # 0 <= RA <= 360
        dec = 30.0   # -90 <= Dec <= 90
        rotSkyPos = 0.0
        self.camera = ComCam()
        self.camera.setObsMetaData(ra, dec, rotSkyPos, mjd=59580.0)

        self.stars = StarData([123, 456, 789], [0.1, 0.2, 0.3],
                              [2.1, 2.2, 2.3], [2.0, 3.0, 4.0],
                              [2.1, 2.1, 4.1], [2.2, 3.2, 4.2],
                              [2.3, 3.3, 4.3], [], [])

    def testCamera(self):

        camera = self.camera
        stars = self.stars

        # Test to get the camera sensor
        # detector = camera.getSensor("center")
        # self.assertEqual(list(detector), ["R:2,2 S:1,1"]) 

        detector = self.camera.getWavefrontSensor()

        # Test to get four camera corners
        corners = detector["R:2,2 S:1,1"]
        self.assertEqual(len(corners), 4)
        
        # Test to transform the stars coordinate to pixel
        stars.setDetector("R:2,2 S:1,1")
        stars = camera.populatePixelFromRADecl(stars)

        self.assertEqual(len(stars.getRaInPixel()), 3)

        # Test to remove stars not on detector
        stars = camera.removeStarsNotOnDetector(stars, 1e7)
        self.assertEqual(len(stars.getRA()), 3)
        self.assertEqual(len(stars.getId()), 3)
        self.assertEqual(len(stars.getRaInPixel()), 3)
        self.assertEqual(len(stars.getDeclInPixel()), 3)
        self.assertEqual(len(stars.getMag(FilterType.U)), 3)
        self.assertEqual(len(stars.getMag(FilterType.G)), 3)
        self.assertEqual(len(stars.getMag(FilterType.R)), 3)
        self.assertEqual(len(stars.getMag(FilterType.I)), 3)
        self.assertEqual(len(stars.getMag(FilterType.Z)), 0)
        self.assertEqual(len(stars.getMag(FilterType.Y)), 0)

        stars = camera.removeStarsNotOnDetector(stars, 0)
        self.assertEqual(stars.getRA().tolist(), [])
        self.assertEqual(stars.getId().tolist(), [])
        self.assertEqual(stars.getRaInPixel().tolist(), [])
        self.assertEqual(stars.getDeclInPixel().tolist(), [])
        self.assertEqual(stars.getMag(FilterType.U).tolist(), [])
        self.assertEqual(stars.getMag(FilterType.G).tolist(), [])
        self.assertEqual(stars.getMag(FilterType.R).tolist(), [])
        self.assertEqual(stars.getMag(FilterType.I).tolist(), [])
        self.assertEqual(stars.getMag(FilterType.Z).tolist(), [])
        self.assertEqual(stars.getMag(FilterType.Y).tolist(), [])

        # Test to get the correct sensor type
        # self.assertEqual(len(camera.getSciCcdList()), 189)
        # self.assertEqual(len(camera.getWfsCcdList()), 8)

        # Test to get the CCD dimension
        self.assertEqual(camera.getCcdDim("R:2,2 S:1,1"), (4072, 4000))


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
