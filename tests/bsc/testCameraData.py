import unittest

from lsst.sims.utils import ObservationMetaData

from lsst.ts.wep.bsc.StarData import StarData
from lsst.ts.wep.bsc.CameraData import ComCam


class TestCameraData(unittest.TestCase):
    """Test the ComCam class."""

    # Boresight (unit: degree)
    RA = 0.0    # 0 <= RA <= 360
    Dec = 30.0   # -90 <= Dec <= 90

    # Camera rotation
    cameraRotation = 0.0
    cameraMJD = 59580.0

    # Get corners of wavefront sensors for this observation field
    obs = ObservationMetaData(pointingRA = RA, pointingDec = Dec, 
                              rotSkyPos = cameraRotation, mjd = cameraMJD)
    camera = None

    # Stars
    stars = None

    def setUp(self):
        self.camera = ComCam()
        self.camera.initializeDetectors()
        self.stars = StarData([123, 456, 789], [0.1, 0.2, 0.3],
                              [2.1, 2.2, 2.3], [2.0, 3.0, 4.0],
                              [2.1, 2.1, 4.1], [2.2, 3.2, 4.2],
                              [2.3, 3.3, 4.3], [2.4, 3.4, 4.4],
                              [2.5, 3.5, 4.5])

    def testCamera(self):

        camera = self.camera
        stars = self.stars

        # Test to get the camera sensor
        detector = camera.getSensor(self.obs, "center")
        self.assertEqual(list(detector), ["R:2,2 S:1,1"]) 

        # Test to get four camera corners
        corners = detector["R:2,2 S:1,1"]
        self.assertEqual(len(corners), 4)
        
        # Test to transform the stars coordinate to pixel
        stars.populateDetector("R:2,2 S:1,1")
        camera.populatePixelFromRADecl(stars, self.obs)

        self.assertEqual(len(stars.getRaInPixel()), 3)

        # Test to remove stars not on detector
        camera.removeStarsNotOnDetectorSimple(stars, self.obs, 1e7)
        self.assertEqual(len(stars.getRA()), 3)
        camera.removeStarsNotOnDetectorSimple(stars, self.obs, 0)
        self.assertEqual(stars.getRA().tolist(), [])

        # Test to get the correct sensor type
        self.assertEqual(len(camera.getSciCcdList()), 189)
        self.assertEqual(len(camera.getWfsCCdList()), 8)

        # Test to get the CCD dimension
        self.assertEqual(camera.getCcdDim("R:2,2 S:1,1"), (4072, 4000))


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
