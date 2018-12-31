import os
import numpy as np
from astropy.io.fits import getheader

from lsst.ts.wep.bsc.Filter import Filter
from lsst.ts.wep.bsc.CamFactory import CamFactory
from lsst.ts.wep.bsc.DatabaseFactory import DatabaseFactory

from lsst.ts.wep.LocalDatabaseDecorator import LocalDatabaseDecorator


class SourceSelector(object):

    UWdb = "UWdb"
    LocalDb = "LocalDb"

    CAMERA_MJD = 59580.0
    STAR_RADIUS_IN_PIXEL = 63
    SPACING_COEFF = 2.5

    def __init__(self, camType, bscDbType):
        """Initialize the source selector class.

        Parameters
        ----------
        camType : CamType
            Camera type.
        bscDbType : BscDbType
            Bright star catalog (BSC) database type.
        """

        self.camera = CamFactory.createCam(camType)
        self.db = DatabaseFactory.createDb(bscDbType)
        self.filter = Filter()

        self.maxDistance = 0.0
        self.maxNeighboringStar = 0

        # Configurate the criteria of neighboring stars
        self.configNbrCriteria(self.STAR_RADIUS_IN_PIXEL, self.SPACING_COEFF,
                               maxNeighboringStar=0)

    def configNbrCriteria(self, starRadiusInPixel, spacingCoefficient,
                          maxNeighboringStar=0):
        """Set the neighboring star criteria to decide the scientific target.

        Parameters
        ----------
        starRadiusInPixel : float, optional
            Diameter of star. For the defocus = 1.5 mm, the star's radius is
            63 pixel.
        spacingCoefficient : float, optional
            Maximum distance in units of radius one donut must be considered
            as a neighbor.
        maxNeighboringStar : int, optional
            Maximum number of neighboring stars. (the default is 0.)
        """

        self.maxDistance = starRadiusInPixel * spacingCoefficient
        self.maxNeighboringStar = int(maxNeighboringStar)

    def connect(self, *kwargs):
        """Connect the database.

        Parameters
        ----------
        *kwargs : str or *list
            Information to connect to the database.
        """

        self.db.connect(*kwargs)

    def disconnect(self):
        """Disconnect the database."""

        self.db.disconnect()

    def setFilter(self, filterType):
        """Set the filter type.

        Parameters
        ----------
        filterType : FilterType
            Filter type.
        """

        self.filter.setFilter(filterType)

    def getFilter(self):
        """Get the filter type.

        Returns
        -------
        FilterType
            Filter type.
        """

        return self.filter.getFilter()

    def setObsMetaData(self, ra, dec, rotSkyPos):
        """Set the observation meta data.

        Parameters
        ----------
        ra : float
            Pointing ra in degree.
        dec : float
            Pointing decl in degree.
        rotSkyPos : float
            The orientation of the telescope in degrees.
        """

        self.camera.setObsMetaData(ra, dec, rotSkyPos, mjd=self.CAMERA_MJD)

    def getTargetStar(self, offset=0):
        """Get the target stars by querying the database.

        Parameters
        ----------
        offset : float, optional
            Offset to the dimension of camera. If the detector dimension is 10
            (assume 1-D), the star's position between -offset and 10+offset
            will be seem to be on the detector. (the default is 0.)

        Returns
        -------
        dict
            Information of neighboring stars and candidate stars with the name
            of sensor as a dictionary.
        dict
            Information of stars with the name of sensor as a dictionary.
        dict
            (ra, dec) of four corners of each sensor with the name
            of sensor as a list. The dictionary key is the sensor name.
        """

        cameraFilter = self.getFilter()
        wavefrontSensors = self.camera.getWavefrontSensor()
        lowMagnitude, highMagnitude = self.filter.getMagBoundary()

        # Query the star database
        starMap = dict()
        neighborStarMap = dict()
        for detector, wavefrontSensor in wavefrontSensors.items():

            # Get stars in this wavefront sensor for this observation field
            stars = self.db.query(cameraFilter, wavefrontSensor[0],
                                  wavefrontSensor[1], wavefrontSensor[2],
                                  wavefrontSensor[3])

            # Set the detector information for the stars
            stars.setDetector(detector)

            # Populate pixel information for stars
            populatedStar = self.camera.populatePixelFromRADecl(stars)

            # Get the stars that are on the detector
            starsOnDet = self.camera.getStarsOnDetector(populatedStar, offset)
            starMap[detector] = starsOnDet

            # Check the candidate of bright stars based on the magnitude
            indexCandidate = starsOnDet.checkCandidateStars(
                                cameraFilter, lowMagnitude, highMagnitude)

            # Determine the neighboring stars based on the distance and
            # allowed number of neighboring stars
            neighborStar = starsOnDet.getNeighboringStar(
                                indexCandidate, self.maxDistance, cameraFilter,
                                self.maxNeighboringStar)
            neighborStarMap[detector] = neighborStar

        return neighborStarMap, starMap, wavefrontSensors

    def insertToBSC(self, neighborStarMap):
        """
        
        Insert the neighboring star data into the local database.
        
        Arguments:
            neighborStarMap {[NeighboringStar]} -- Information of neighboring stars.
        
        Raises:
            ValueError -- Not the local database.
        """
        
        # Check the database is the local database or not
        if (self.name != self.LocalDb):
            raise ValueError("Can not insert data into '%s'." % self.name)

        # Insert the data 
        for detector, singleNeighborStarMap in neighborStarMap.items():
            self.db.insertData(self.getFilter(), singleNeighborStarMap)

    def generateBSC(self, localDb):
        """
        
        Generate the bright star catalog.
        
        Arguments:
            localDb {[database]} -- Local database to put the bright star catalog.
        
        Raises:
            ValueError -- Not remote UW database.
            TypeError -- Not ComCam type camera.
        """

        # Check the database is the UW database or not
        if (self.name != self.UWdb):
            raise ValueError("Can not generate BSC from '%s'." % self.name)

        # Check the camera is comcam or not
        if (not isinstance(self.camera, ComCam)):
            raise TypeError("Camera should be ComCam type.")

        # Boresight (unit: degree)
        delta = 0.2
        RaArray = np.arange(0, 360, delta)
        DecArray = np.append(np.arange(-90, 90, delta), 90)

        # Set the rotation angle
        cameraRotation = 0.0

        # Set the filter in localDb to be the same as the remote UW database
        localDb.setFilter(self.getFilter())

        # Go through all (RA, Dec)
        for RA in RaArray:
            for Dec in DecArray:
                # Do the query
                neighborStarMap, starMap, wavefrontSensors = self.getTargetStar((RA, Dec), cameraRotation, 
                                                                                orientation="center", 
                                                                                offset=self.maxDistance)

                # Write data into the local database
                localDb.insertToBSC(neighborStarMap)

    def searchRaDecl(self, ra, decl):
        """
        
        Search the star id based on ra, decl.
        
        Arguments:
            ra {[float]} -- ra in degree (0 deg - 360 deg).
            decl {[float]} -- decl in degree (-90 deg - 90 deg).
        
        Returns:
            [int] -- Star ID in database. It is noted that the id will be "simobjid" if the database is 
                     the remote UW database.
        """

        starID = []

        # Get the star id from the database
        if (self.name == self.UWdb):
            starID = self.db.searchRaDecl(self.tableName, ra, decl)

            # Change the data type from decimal to int to keep the same data type as local database
            # This might be removed when the local database switchs to mssql.
            starID.append((int(starID.pop()[0]),))

        elif (self.name == self.LocalDb):
            starID = self.db.searchRaDecl(self.filter.getFilter(), ra, decl)
            
        return starID

    def updateBSC(self, listID, listOfItemToChange, listOfNewValue):
        """
        
        Update data based on the id.
        
        Arguments:
            listID {[int]} -- ID list to change.
            listOfItemToChange {[string]} -- Item list (simobjid, ra, decl, mag, bright_star) to change.
            listOfNewValue {[valueType]} -- New value list.
        
        Raises:
            ValueError -- Not local database.
        """

        # Check the database is the local database or not
        if (self.name != self.LocalDb):
            raise ValueError("Can not update BSC in '%s'." % self.name)

        # Update the bright star catalog
        self.db.updateData(self.filter.getFilter(), listID, listOfItemToChange, listOfNewValue)

    def trimMargin(self, neighborStarMap, trimMarginInPixel):
        """
        
        Trim the candidate stars if they or related neighboring stars are outside of boundary. 
        
        Arguments:
            neighborStarMap{[list]} -- Information of neighboring stars and candidate stars with 
                                       the name of sensor as a list.
            trimMarginInPixel {[float]} -- Trimed boundary in pixel. if the ccd dimension is (d1, d2), only stars 
                                           inside (trimMarginInPixel < x1 < d1-trimMarginInPixel) and 
                                           (trimMarginInPixel < x2 < d2-trimMarginInPixel) will be left.
        
        Raises:
            ValueError -- trimMarginInPixel is less than 0.
            ValueError -- trimMarginInPixel is bigger than the half of CCD dimension.
        """

        # Check the boundary of trimMarginInPixel
        if (trimMarginInPixel<0):
            raise ValueError("The trimmed boudary pixel < 0.")

        # Trim the stars that are in the margin.
        for detector, neighborStar in neighborStarMap.items():

            trimmedCandidateStarNum = 0

            # Detector dimension
            dim1, dim2 = self.camera.getCcdDim(detector)

            # Check the boundary of trimMarginInPixel
            if (trimMarginInPixel > min(dim1, dim2)/2):
                raise ValueError("trimMarginInPixel ('%f') >= half of CCD's dimension." % trimMarginInPixel)

            # Copy a new dictionary to avoid the iteration error for changing size of iteration
            neighborStarSimobjID = neighborStar.SimobjID.copy()

            # Use the candidate star as the unit to check stars are inside the boundary or not
            for candidateStar, neighboringStars in neighborStarSimobjID.items():

                # Get all stars (candidateStar: string + neighboringStars: list) in this item
                # Use the List[:] to get the copy of list
                allStars = neighboringStars[:]
                allStars.append(candidateStar)

                needToTrim = False
                # Check the coordinate for each star
                for star in allStars:
                    coord1, coord2 = neighborStar.RaDeclInPixel[star]

                    # Check the star inside the boundary or not
                    if (coord1 <= trimMarginInPixel or coord1 >= dim1-trimMarginInPixel or 
                        coord2 <= trimMarginInPixel or coord2 >= dim2-trimMarginInPixel):
                        needToTrim = True
                        break
                
                # The candidate/ neighboring stars are outside of tbe boundary
                if (needToTrim):
                    trimmedCandidateStarNum += 1

                    # Add "None" to avoid the raised error that there is the unfound key
                    neighborStar.SimobjID.pop(candidateStar, None)

            if (trimmedCandidateStarNum != 0):
                print("Trimmed candidate stars on %s: %d." % (detector, trimmedCandidateStarNum))

def calcPixPos(fitsFilePath, raList, decList, extLayer=0):
    """
    
    Calculate the pixel positions based on the FITS header data. The working formula
    is provided by John.
    
    Arguments:
        fitsFilePath {[str]} -- FITS file path.
        raList {[list]} -- List of RA in degree.
        decList {[list]} -- List of Dec in degree.
    
    Keyword Arguments:
        extLayer {int} -- Extension layer in degree. (default: {0})
    
    Returns:
        [list] -- List of x position in pixel.
        [list] -- List of y position in pixel.
    """

    # Get the header data
    hdr = getheader(fitsFilePath, int(extLayer))

    # Get the needed parameters
    crval1 = hdr["CRVAL1"]
    crval2 = hdr["CRVAL2"]
    cd1_1 = hdr["CD1_1"]
    cd1_2 = hdr["CD1_2"]
    cd2_1 = hdr["CD2_1"]
    cd2_2 = hdr["CD2_2"]
    crpix1 = hdr["CRPIX1"]
    crpix2 = hdr["CRPIX2"]

    # Calculate the pixel position
    xPosList = []
    yPosList = []
    for ra, dec in zip(raList, decList):

        # Change the direction. RA is in (-180, 180).
        if (ra > 180):
            ra = ra - 360

        # Calculate the pixel position based on the formula provided by John
        xPos = crpix1 + ( (ra-crval1)*cd2_2 - (dec-crval2)*cd1_2 ) / (cd2_2*cd1_1 - cd1_2*cd2_1)
        yPos = crpix2 + ( (ra-crval1)*cd2_1 - (dec-crval2)*cd1_1 ) / (cd1_2*cd2_1 - cd2_2*cd1_1)

        xPosList.append(xPos)
        yPosList.append(yPos)

    return xPosList, yPosList


if __name__ == "__main__":
    pass
