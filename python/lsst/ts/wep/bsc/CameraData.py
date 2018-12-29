import numpy as np


from lsst.sims.coordUtils.CameraUtils import raDecFromPixelCoords, \
                                             pixelCoordsFromRaDec

from lsst.sims.utils import ObservationMetaData

from lsst.ts.wep.bsc.WcsSol import WcsSol

from lsst.afw.cameraGeom import WAVEFRONT, SCIENCE

from lsst.ts.wep.bsc.StarData import StarData

from lsst.ts.wep.Utility import FilterType

# Ignore the warning of "WARNING: ErfaWarning: ERFA function "taiutc" yielded 1 
# of "dubious year (Note 4)" [astropy._erfa.core]"

class CameraData(object):
    
    def __init__(self, cameraCollection):
        """
        
        Initiate the camera for bright star catalog to use.
        
        Arguments:
            cameraCollection {[camera]} -- A collection of detectors that also supports 
                                           coordinate transformation.
        """

        self._camera = cameraCollection

        # Dictionary of (x, y) coordinates of detector corners and dimensions of detection
        self._corners = {}
        self._dimension = {}

        # List of camera CCD
        self._wfsCcd = []
        self._sciCcd = []

    def getCameraCollection(self):
        """
        
        Get the camera collection.
        
        Returns:
            [camera] -- A collection of detectors that also supports coordinate transformation.
        """
        
        return self._camera
        
    def initializeDetectors(self):
        """
        Initializes the camera wavefront detectors.
        """

        for detector in self._camera:
                if detector.getType() in (WAVEFRONT, SCIENCE):

                    detectorName = detector.getName()

                    # Collect the ccd name
                    if (detector.getType() == WAVEFRONT):
                        self._wfsCcd.append(detectorName)
                    elif (detector.getType() == SCIENCE):
                        self._sciCcd.append(detectorName)

                    bbox = detector.getBBox()
                    xmin = bbox.getMinX()
                    xmax = bbox.getMaxX()
                    ymin = bbox.getMinY()
                    ymax = bbox.getMaxY()
                    self._corners[detectorName] = (np.array([xmin, xmin, xmax, xmax]), 
                                                    np.array([ymin, ymax, ymin, ymax]))

                    # The CCD dimension here is an estimation. 
                    # Based on LCA-13381, there are three types of sensors.
                    # e2V CCD250: 40.04 mm x 40.96 mm
                    # STA 4400: 20.00 mm x 40.72 mm
                    # STA 3800C: 40.00 mm x 40.72 mm

                    dim1, dim2 = bbox.getDimensions()
                    self._dimension[detectorName] = (int(dim1), int(dim2))  

    def populatePixelFromRADecl(self, stars, obs):
        """
        
        Populates the RAInPixel and DeclInPixel coordinates in the StarData stars using the lsst-sims 
        stack.
        
        Arguments:
            stars {[StarData]} -- The stars to populate.
            obs {[metadata]} -- The observation meta data (found in the lsst-sims stack) that defines 
                                the pointing.
        """

        ra = stars.getRA()
        decl = stars.getDecl()
        raInPixel, declInPixel = pixelCoordsFromRaDec(ra = ra, dec = decl, obs_metadata = obs,
                                                      epoch = 2000.0, 
                                                      chipName = np.array([stars.getDetector()] * len(ra)), 
                                                      camera = self._camera, includeDistortion = True)
        stars.setRaInPixel(raInPixel)
        stars.setDeclInPixel(declInPixel)
        
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
        
        keep = [index for index in range(len(stars.getRA())) 
                if stars.getRaInPixel()[index] >= -offset and stars.getRaInPixel()[index] <= self._dimension[stars.getDetector()][0]+offset 
                and stars.getDeclInPixel()[index] >= -offset and stars.getDeclInPixel()[index] <= self._dimension[stars.getDetector()][1]+offset]
        
        starsRA = [stars.getRA()[index] for index in keep]
        stars.setRA(starsRA)

        starsRAInPixel = [stars.getRaInPixel()[index] for index in keep]
        stars.setRaInPixel(starsRAInPixel)
     
        starsDecl = [stars.getDecl()[index] for index in keep]
        stars.setDecl(starsDecl)

        stars.DeclInPixel = [stars.getDeclInPixel()[index] for index in keep]
        stars.setDeclInPixel(starsDecl)
        
        # Check the empty information
        if (len(stars.getMag(FilterType.U)) != 0):
            starsLSSTMagU = [stars.getMag(FilterType.U)[index] for index in keep]
            stars.setMag(FilterType.U, starsLSSTMagU)
         
        if (len(stars.getMag(FilterType.G)) != 0):
            starsLSSTMagG = [stars.getMag(FilterType.G)[index] for index in keep]
            stars.setMag(FilterType.G, starsLSSTMagG)

        if (len(stars.getMag(FilterType.R)) != 0):
            starsLSSTMagR = [stars.getMag(FilterType.R)[index] for index in keep]
            stars.setMag(FilterType.R, starsLSSTMagR)

        if (len(stars.getMag(FilterType.I)) != 0):
            starsLSSTMagI = [stars.getMag(FilterType.I)[index] for index in keep]
            stars.setMag(FilterType.I, starsLSSTMagI)

        if (len(stars.getMag(FilterType.Z)) != 0):
            starsLSSTMagZ = [stars.getMag(FilterType.Z)[index] for index in keep]
            stars.setMag(FilterType.Z, starsLSSTMagZ)

        if (len(stars.getMag(FilterType.Y)) != 0):
            starsLSSTMagY = [stars.getMag(FilterType.Y)[index] for index in keep]
            stars.setMag(FilterType.Y, starsLSSTMagY)

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

            coords = self._corners[detector]

            ra, dec = raDecFromPixelCoords(coords[0], coords[1], [detector]*len(coords[0]),
                                           camera=self._camera, obs_metadata=obs,
                                           epoch=2000.0, includeDistortion=True)   

            ra_dec_out[detector] = [(ra[0], dec[0]), (ra[1], dec[1]), (ra[2], dec[2]), (ra[3], dec[3])]

        return ra_dec_out

    def getWfsCCdList(self):
        """
        
        Get the list of wavefront sensor list.
        
        Returns:
            [list] -- CCD list.
        """

        return self._wfsCcd

    def getSciCcdList(self):
        """
        
        Get the list of science sensor list.
        
        Returns:
            [list] -- CCD list.
        """

        return self._sciCcd

    def getCcdDim(self, detectorName):
        """
        
        Get the CCD dimension.
        
        Arguments:
            detectorName {[string]} -- Detector Name.
        
        Returns:
            [tuple] -- CCD dimension in pixel.
        """

        return self._dimension[detectorName]

    def getWavefrontSensor(self):
        raise NotImplementedError("Subclass must implement the abstract method.")


if __name__ == "__main__":
    pass
