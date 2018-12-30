import numpy as np

from lsst.ts.wep.bsc.WcsSol import WcsSol
from lsst.ts.wep.Utility import FilterType


class CameraData(object):
    
    def __init__(self, camera):
        """Initialize the camera data class.

        Parameters
        ----------
        camera : lsst.afw.cameraGeom.camera.camera.Camera
            A collection of Detectors that also supports coordinate
            transformation. (the default is None.)
        """

        self._wcs = WcsSol(camera=camera)

        # List of wavefront sensor CCD name
        self._wfsCcd = []

        # Dictionary of (x, y) coordinates of detector corners and dimensions
        # of detection. The dictonary key is the ccd name.
        self._corners = dict()
        self._dimension = dict()

    def setWfsCcdList(self, wfsCcdList):

        self._wfsCcd = wfsCcdList

    def setWfsCorners(self, wfsCorners):

        self._corners = wfsCorners

    def getWfsCorners(self):

        return self._corners

    def setCcdDims(self, ccdDims):

        self._dimension = ccdDims

    def _initDetectors(self, detectorType):
        """Initializes the camera detectors.

        Parameters
        ----------
        detectorType : lsst.afw.cameraGeom.detector.detector.DetectorType
            Detector type.
        """

        for detector in self._wcs.getCamera():

            if (detector.getType() == detectorType):

                # Collect the ccd name
                detectorName = detector.getName()
                self._wfsCcd.append(detectorName)

                # Get the detector corners
                bbox = detector.getBBox()
                xmin = bbox.getMinX()
                xmax = bbox.getMaxX()
                ymin = bbox.getMinY()
                ymax = bbox.getMaxY()
                self._corners[detectorName] = \
                                (np.array([xmin, xmin, xmax, xmax]), 
                                 np.array([ymin, ymax, ymin, ymax]))

                # The CCD dimension here is an estimation. 
                # Based on LCA-13381, there are three types of sensors.
                # e2V CCD250: 40.04 mm x 40.96 mm
                # STA 4400: 20.00 mm x 40.72 mm
                # STA 3800C: 40.00 mm x 40.72 mm

                dim1, dim2 = bbox.getDimensions()
                self._dimension[detectorName] = (int(dim1), int(dim2))

    def setObsMetaData(self, ra, dec, rotSkyPos, mjd=59580.0):
        """Set the observation meta data.

        Parameters
        ----------
        ra : float
            Pointing ra in degree.
        dec : float
            Pointing decl in degree.
        rotSkyPos : float
            The orientation of the telescope in degrees.
        mjd : float
            Camera MJD. (the default is 59580.0.)
        """

        self._wcs.setObsMetaData(ra, dec, rotSkyPos, mjd=mjd)

    def populatePixelFromRADecl(self, stars):
        """Populates the RAInPixel and DeclInPixel coordinates to the stars.

        Parameters
        ----------
        stars : StarData
            The stars to populate.

        Returns
        -------
        StarData
            The stars with x-, y-pixel data populated.
        """

        ra = stars.getRA()
        decl = stars.getDecl()
        chipName = np.array([stars.getDetector()] * len(ra))
        raInPixel, declInPixel = self._wcs.pixelCoordsFromRaDec(
            ra, decl, chipName=chipName, epoch=2000.0, includeDistortion=True)

        stars.setRaInPixel(raInPixel)
        stars.setDeclInPixel(declInPixel)

        return stars

    def removeStarsNotOnDetector(self, stars, offset):
        """Remove the 'stars' from the input stars that are not on the detector
        according to the pixel position.

        Parameters
        ----------
        stars : StarData
            Star information.
        offset : float
            The offset to dimension of camera. This is for generating the local
            database of bright star catalog for the condition that the bright
            star is near the edge of ccd.

        Returns
        -------
        StarData
            The stars object without the removed 'stars' that are not on the
            detector.
        """

        # Get the index that will keep the data
        starsRaInPixel = stars.getRaInPixel()
        starsDeclInPixel = stars.getDeclInPixel()
        ccdDim = self.getCcdDim(stars.getDetector())

        keep = []
        for ii in range(len(starsRaInPixel)):
            if (-offset <= starsRaInPixel[ii] <= ccdDim[0] + offset and 
                -offset <= starsDeclInPixel[ii] <= ccdDim[1] + offset):
                keep.append(ii)

        # Remove the stars that are not on the detector
        stars.setId(self._getKeepItem(stars.getId(), keep))
        stars.setRA(self._getKeepItem(stars.getRA(), keep))
        stars.setDecl(self._getKeepItem(stars.getDecl(), keep))
        stars.setRaInPixel(self._getKeepItem(starsRaInPixel, keep))
        stars.setDeclInPixel(self._getKeepItem(starsDeclInPixel, keep))

        for filterType in FilterType:
            magArray = stars.getMag(filterType)
            if (len(magArray) != 0):
                stars.setMag(filterType, self._getKeepItem(magArray, keep))

        return stars

    def _getKeepItem(self, valArray, keep):
        """Get the keep items in array.

        Parameters
        ----------
        valArray : list or 1-D numpy.ndarray
            Value array.
        keep : list[int]
            List of keep index.

        Returns
        -------
        list
            The array that has only the keep values.
        """

        return [valArray[idx] for idx in keep]

    def getWfsCcdList(self):
        """Get the list of wavefront sensor CCD list.

        Returns
        -------
        list
            Wavefront sensor CCD list.
        """

        return self._wfsCcd

    def getWavefrontSensor(self):
        """

        Get the corners of LSST curvature wavefront sensors in (ra, dec) based
        on the wavefront sensor list.

        Returns:
            [dict] -- (ra, dec) of four corners of each sensor with the name
            of sensor as a list. The dictionary key is the sensor name.
        """

        return self._getDetectorRaDec(self._wfsCcd)

    def _getDetectorRaDec(self, detectorList):
        """Get the (ra, dec) of CCD corners in the detector list.

        Parameters
        ----------
        detectorList : list
            List of detectors. For example, ["R:2,2 S:1,1", "R:2,2 S:0,1"].

        Returns
        -------
        dict
            This method returns a dict of list.  The dict is keyed on the name
            of the wavefront sensor.  The list contains the (RA, Dec)
            coordinates of the corners of that sensor (RA, Dec are paired as
            tuples). For example, output['R:0,0 S:2,2B'] = [(23.0, -5.0),
            (23.1, -5.0), (23.0, -5.1), (23.1, -5.1)] would mean that the
            wavefront sensor named 'R:0,0 S:2,2B' has its corners at RA 23,
            Dec -5; RA 23.1, Dec -5; RA 23, Dec -5.1; and RA 23.1, Dec -5.1 
            Coordinates are in degrees.
        """

        ra_dec_out = dict()
        for detector in detectorList:

            coords = self._corners[detector]
            xPix = coords[0]
            yPix = coords[1]

            chipName = np.array([detector] * len(xPix))
            ra, dec = self._wcs.raDecFromPixelCoords(
                xPix, yPix, chipName, epoch=2000.0, includeDistortion=True)

            ra_dec_out[detector] = [(ra[0], dec[0]), (ra[1], dec[1]),
                                    (ra[2], dec[2]), (ra[3], dec[3])]

        return ra_dec_out

    def getCcdDim(self, detectorName):
        """Get the CCD dimension.

        Parameters
        ----------
        detectorName : str
            Detector Name (e.g. "R:2,2 S:1,1").

        Returns
        -------
        tuple
            CCD dimension in pixel.
        """

        return self._dimension[detectorName]


if __name__ == "__main__":
    pass
