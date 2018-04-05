import numpy as np
import unittest

# instantiate the LSST camera model
from lsst.obs.lsstSim import LsstSimMapper

from lsst.sims.coordUtils.CameraUtils import raDecFromPixelCoords, pixelCoordsFromRaDec
from lsst.sims.utils import ObservationMetaData
from lsst.afw.cameraGeom import WAVEFRONT, SCIENCE

from lsst.ts.wep.bsc.StarData import StarData

# Ignore the warning of "WARNING: ErfaWarning: ERFA function "taiutc" yielded 1 
# of "dubious year (Note 4)" [astropy._erfa.core]"

class CameraData(object):

    LSST = "lsst"
    COMCAM = "comcam"
    
    def __init__(self, name, cameraCollection):
        """
        
        Initiate the camera for bright star catalog to use.
        
        Arguments:
            name {[string]} -- Name of camera.
            cameraCollection {[camera]} -- A collection of detectors that also supports 
                                           coordinate transformation.
        """

        self.name = name
        self.__camera = cameraCollection

        # Dictionary of (x, y) coordinates of detector corners and dimensions of detection
        self.__corners = {}
        self.__dimension = {}

        # List of camera CCD
        self.__wfsCcd = []
        self.__sciCcd = []

    def getCameraCollection(self):
        """
        
        Get the camera collection.
        
        Returns:
            [camera] -- A collection of detectors that also supports coordinate transformation.
        """
        
        return self.__camera
        
    def initializeDetectors(self):
        """
        Initializes the camera wavefront detectors.
        """

        for detector in self.__camera:
                if detector.getType() in (WAVEFRONT, SCIENCE):

                    detectorName = detector.getName()

                    # Collect the ccd name
                    if (detector.getType() == WAVEFRONT):
                        self.__wfsCcd.append(detectorName)
                    elif (detector.getType() == SCIENCE):
                        self.__sciCcd.append(detectorName)

                    bbox = detector.getBBox()
                    xmin = bbox.getMinX()
                    xmax = bbox.getMaxX()
                    ymin = bbox.getMinY()
                    ymax = bbox.getMaxY()
                    self.__corners[detectorName] = (np.array([xmin, xmin, xmax, xmax]), 
                                                    np.array([ymin, ymax, ymin, ymax]))

                    # The CCD dimension here is an estimation. 
                    # Based on LCA-13381, there are three types of sensors.
                    # e2V CCD250: 40.04 mm x 40.96 mm
                    # STA 4400: 20.00 mm x 40.72 mm
                    # STA 3800C: 40.00 mm x 40.72 mm

                    dim1, dim2 = bbox.getDimensions()
                    self.__dimension[detectorName] = (int(dim1), int(dim2))  

    def populatePixelFromRADecl(self, stars, obs):
        """
        
        Populates the RAInPixel and DeclInPixel coordinates in the StarData stars using the lsst-sims 
        stack.
        
        Arguments:
            stars {[StarData]} -- The stars to populate.
            obs {[metadata]} -- The observation meta data (found in the lsst-sims stack) that defines 
                                the pointing.
        """

        ra = stars.RA
        decl = stars.Decl
        raInPixel, declInPixel = pixelCoordsFromRaDec(ra = ra, dec = decl, obs_metadata = obs,
                                                      epoch = 2000.0, 
                                                      chipName = np.array([stars.Detector] * len(stars.RA)), 
                                                      camera = self.__camera, includeDistortion = True)
        stars.populateRAData(raInPixel)
        stars.populateDeclData(declInPixel)
        
    def removeStarsNotOnDetectorSimple(self, stars, obs, offset):
        """
        
        Removes the stars from the StarData stars that are not on the detector using pixel data.

        Arguments:
            stars {[StarData]} -- Star information.
            obs {[metadata]} -- The observation meta data (found in the lsst-sims stack) that defines 
                                the pointing.
            offset {[float]} -- The offset to dimension of camera. This is for generating the local 
                                database of bright star catalog for the condition that the bright star 
                                is near the edge of ccd.
        """
        
        keep = [index for index in range(len(stars.RA)) 
                if stars.RAInPixel[index] >= -offset and stars.RAInPixel[index] <= self.__dimension[stars.Detector][0]+offset 
                and stars.DeclInPixel[index] >= -offset and stars.DeclInPixel[index] <= self.__dimension[stars.Detector][1]+offset]
        
        stars.RA = [stars.RA[index] for index in keep]
        stars.RAInPixel = [stars.RAInPixel[index] for index in keep]
        stars.Decl = [stars.Decl[index] for index in keep]
        stars.DeclInPixel = [stars.DeclInPixel[index] for index in keep]
        
        # Check the empty information
        if (stars.LSSTMagU):
            stars.LSSTMagU = [stars.LSSTMagU[index] for index in keep]
         
        if (stars.LSSTMagG):
            stars.LSSTMagG = [stars.LSSTMagG[index] for index in keep]
         
        if (stars.LSSTMagR):
            stars.LSSTMagR = [stars.LSSTMagR[index] for index in keep]
         
        if (stars.LSSTMagI):
            stars.LSSTMagI = [stars.LSSTMagI[index] for index in keep]
         
        if (stars.LSSTMagZ):
            stars.LSSTMagZ = [stars.LSSTMagZ[index] for index in keep]
         
        if (stars.LSSTMagY):
            stars.LSSTMagY = [stars.LSSTMagY[index] for index in keep]
    
    def getDetectorRaDec(self, camera_mapper, obs):
        """
        
        Get the (ra, dec) of ccd corners.

        Arguments:
            camera_mapper {[metadata]} -- camera_mapper is the sensor ID on LSST camera map. 
                                          The detail of map is at:
                    https://confluence.lsstcorp.org/display/LSWUG/Representation+of+a+Camera
                                          The format looks like 'R:2,2 S:2,0' for science sensor 
                                          and 'R:0,0 S:2,2B' for corner wavefront sensor. 
            obs {[metadata]} -- Instantiation of ObservationMetaData that describes the pointing
                                of the telescope.
        
        Returns:
            [list] -- This method returns a dict of list.  The dict is keyed on the name of the
                      wavefront sensor.  The list contains the (RA, Dec) coordinates of the corners
                      of that sensor (RA, Dec are paired as tuples). For example, 
                      output['R:0,0 S:2,2B'] = [(23.0, -5.0), (23.1, -5.0), (23.0, -5.1), (23.1, -5.1)]
                      would mean that the wavefront sensor named 'R:0,0 S:2,2B' has its corners at
                      RA 23, Dec -5; RA 23.1, Dec -5; RA 23, Dec -5.1; and RA 23.1, Dec -5.1 
                      Coordinates are in degrees.
        """

        ra_dec_out = {}

        for detector in camera_mapper:

            coords = self.__corners[detector]

            ra, dec = raDecFromPixelCoords(coords[0], coords[1], [detector]*len(coords[0]),
                                           camera=self.__camera, obs_metadata=obs,
                                           epoch=2000.0, includeDistortion=True)   

            ra_dec_out[detector] = [(ra[0], dec[0]), (ra[1], dec[1]), (ra[2], dec[2]), (ra[3], dec[3])]

        return ra_dec_out

    def getWfsCCdList(self):
        """
        
        Get the list of wavefront sensor list.
        
        Returns:
            [list] -- CCD list.
        """

        return self.__wfsCcd

    def getSciCcdList(self):
        """
        
        Get the list of science sensor list.
        
        Returns:
            [list] -- CCD list.
        """

        return self.__sciCcd

    def getCcdDim(self, detectorName):
        """
        
        Get the CCD dimension.
        
        Arguments:
            detectorName {[string]} -- Detector Name.
        
        Returns:
            [tuple] -- CCD dimension in pixel.
        """

        return self.__dimension[detectorName]

    def getWavefrontSensor(self):
        raise NotImplementedError("Subclass must implement the abstract method.")

class LsstCamera(CameraData):

    def __init__(self):
        super(LsstCamera, self).__init__(self.LSST, LsstSimMapper().camera)

    def getWavefrontSensor(self, obs):
        """
        
        Get the corners of LSST curvature wavefront sensors in (ra, dec) based on the camera_mapper
        list below.
        
        Arguments:
            obs {[metadata]} -- Instantiation of ObservationMetaData that describes the pointing
                                of the telescope.
        
        Returns:
            [list] -- (ra, dec) of four corners of each sensor with the name of sensor as a list
        """

        # Camera object
        camera_mapper = self.getWfsCCdList()
        ra_dec_out = self.getDetectorRaDec(camera_mapper, obs)

        return ra_dec_out

    def getScineceSensor(self, obs):
        """
        
        Get the corners of LSST science sensors in (ra, dec) based on the camera_mapper list below.
        
        Arguments:
            obs {[metadata]} -- Instantiation of ObservationMetaData that describes the pointing
                                of the telescope.
        
        Returns:
            [list] -- (ra, dec) of four corners of each sensor with the name of sensor as a list
        """

        # Camera object
        camera_mapper = self.getSciCcdList()
        ra_dec_out = self.getDetectorRaDec(camera_mapper, obs)

        return ra_dec_out

class ComCam(CameraData):

    def __init__(self):
        # The comcam's configuration here is approximated by taking the central 
        # raft of lsst camera.    
        super(ComCam, self).__init__(self.COMCAM, LsstSimMapper().camera)

    def getWavefrontSensor(self, obs):
        """
        
        Get the corner sensors of Comcam in (ra, dec) based on the camera_mapper list below.
        The reference is at:
        https://confluence.lsstcorp.org/display/LSWUG/Representation+of+a+Camera
        
        Arguments:
            obs {[metadata]} -- Instantiation of ObservationMetaData that describes the pointing
                                of the telescope.
        
        Returns:
            [list] -- (ra, dec) of four corners of each sensor with the name of sensor as a list
        """

        return self.getSensor(obs, "corner")
    
    def getSensor(self, obs, orientation):
        """
        
        Get the sensors of Comcam in (ra, dec) based on the camera_mapper list below.
        The reference is at:
        https://confluence.lsstcorp.org/display/LSWUG/Representation+of+a+Camera
        
        Arguments:
            obs {[metadata]} -- Instantiation of ObservationMetaData that describes the pointing
        of the telescope.
            orientation {[string]} -- Orientation of camera to decide which sensor to use.
        
        Returns:
            [list] -- (ra, dec) of four corners of each sensor with the name of sensor as a list.
        """

        # Camera object
        if (orientation == "center"):
            camera_mapper = ["R:2,2 S:1,1"]
        elif (orientation == "corner"):
            camera_mapper = ["R:2,2 S:0,2", "R:2,2 S:2,2", "R:2,2 S:0,0", "R:2,2 S:2,0"]
        elif (orientation == "all"):
            camera_mapper = ["R:2,2 S:0,2", "R:2,2 S:1,2", "R:2,2 S:2,2", "R:2,2 S:0,1", 
                             "R:2,2 S:1,1", "R:2,2 S:2,1", "R:2,2 S:0,0", "R:2,2 S:1,0", 
                             "R:2,2 S:2,0"]

        ra_dec_out = self.getDetectorRaDec(camera_mapper, obs)

        return ra_dec_out

class CameraDataTest(unittest.TestCase):
    """
    Test the camera functions.
    """

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
        self.stars = StarData([123, 456, 789], [0.1, 0.2, 0.3], [2.1, 2.2, 2.3], [2.0, 3.0, 4.0], 
                              [2.1, 2.1, 4.1], [2.2, 3.2, 4.2], [2.3, 3.3, 4.3], [2.4, 3.4, 4.4], 
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

        self.assertEqual(len(stars.RAInPixel), 3)

        # Test to remove stars not on detector
        camera.removeStarsNotOnDetectorSimple(stars, self.obs, 1e7)
        self.assertEqual(len(stars.RA), 3)
        camera.removeStarsNotOnDetectorSimple(stars, self.obs, 0)
        self.assertEqual(stars.RA, [])

        # Test to get the correct sensor type
        self.assertEqual(len(camera.getSciCcdList()), 189)
        self.assertEqual(len(camera.getWfsCCdList()), 8)

        # Test to get the CCD dimension
        self.assertEqual(camera.getCcdDim("R:2,2 S:1,1"), (4072, 4000))

if __name__ == "__main__":
 
    # Do the unit test
    unittest.main()



