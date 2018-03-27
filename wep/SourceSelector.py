import os, unittest
import numpy as np
from astropy.io.fits import getheader

from lsst.sims.utils import ObservationMetaData

from bsc.BrightStarDatabase import BrightStarDatabase
from bsc.CameraData import LsstCamera, ComCam 
from bsc.Filter import Filter

from wep.LocalDatabaseDecorator import LocalDatabaseDecorator

class SourceSelector(object):

    UWdb = "UWdb"
    LocalDb = "LocalDb"

    LSST = "lsst"
    COMCAM = "comcam"

    def __init__(self):
        """
        
        Initialize the SourceSelector class.
        """

        self.db = None
        self.tableName = None
        self.name = None

        self.camera = None
        self.cameraMJD = None

        self.maxDistance = np.nan
        self.maxNeighboringStar = np.nan

        self.filter = Filter()

    def setDb(self, dbType):
        """
        
        Set the database.
        
        Arguments:
            dbType {[str]} -- Type of database ("UWdb" or "LocalDb").
        """

        if (dbType == self.UWdb):
            db = BrightStarDatabase()
            tableName = "bright_stars"
        
        elif (dbType == self.LocalDb):
            db = LocalDatabaseDecorator()
            tableName = "BrightStarCatalog"
        else:
            raise ValueError("No '%s' database." % dbType)

        self.db = db
        self.tableName = tableName
        self.name = dbType

    def setCamera(self, cameraType, cameraMJD=59580.0):
        """
        
        Set the camera.
        
        Arguments:
            cameraType {[str]} -- Type of camera ("lsst" or "comcam"). (default: {None})
        
        Keyword Arguments:
            cameraMJD {float} -- Camera MJD. (default: {59580.0})
        """

        # Set the camera and do the initialization
        if (cameraType == self.LSST):
            camera = LsstCamera()
        elif (cameraType == self.COMCAM):
            camera = ComCam()
        else:
            raise ValueError("No '%s' camera." % cameraType)

        self.camera = camera
        self.camera.initializeDetectors()

        # Set the camera MJD
        self.cameraMJD = cameraMJD

    def configSelector(self, cameraType=None, dbType=None, aFilter=None, cameraMJD=59580.0):
        """
        
        Configurate the source selector.
        
        Keyword Arguments:
            cameraType {[str]} -- Type of camera ("lsst" or "comcam"). (default: {None})
            dbType {[str]} -- Type of database ("UWdb" or "LocalDb"). (default: {None})
            aFilter {[str]} -- Active filter type ("u", "g", "r", "i", "z", "y"). 
                                (default: {None})
            cameraMJD {float} -- Camera MJD. (default: {59580.0})
        
        Raises:
            ValueError -- No database type.
            ValueError -- No camera type.
        """

        # Set the data base
        if (dbType is not None):
            self.setDb(dbType)

        # Set the camera mapper
        if (cameraType is not None):
            self.setCamera(cameraType, cameraMJD=cameraMJD)

        # Set the filter
        if (aFilter is not None):
            self.setFilter(aFilter)

    def connect(self, *kwargs):
        """
        
        Connect the database.
        
        Arguments:
            *kwargs {[string]} -- Information to connect to the database.
        """

        self.db.connect(*kwargs)

    def disconnect(self):
        """
        
        Disconnect the database.
        """

        self.db.disconnect()

    def getTargetStar(self, pointing, cameraRotation, orientation=None, offset=0, tableName=None):
        """
        
        Get the target stars by querying the database.
        
        Arguments:
            pointing {[tuple]} -- Camera boresight (RA, Decl) in degree.
            cameraRotation {[float]} -- Camera rotation angle in degree.
        
        Keyword Arguments:
            orientation {[str]} -- Orientation of wavefront sensor(s) on camera. (default: {None})
            offset {[float]} -- Add offset in pixel to sensor dimension for judging stars on detector or not. 
                                offset=0 for normal use. 
                                offset=maxDistance to generate the local database. (default: {0})
            tableName {[str]} -- Table name. (default: {None})
        
        Returns:
            {[dict]} -- Information of neighboring stars and candidate stars with the name of 
                        sensor as a dictionary.
            {[dict]} -- Information of stars with the name of sensor as a dictionary.
            {[dict]} -- Corners of sensor with the name of sensor as a dictionary.
        """

        # Filter type
        cameraFilter = self.filter.getFilter()

        # Regenerate the table name for local database
        if (tableName is None):
            if (self.name == self.LocalDb):
                tableName = self.tableName + self.filter.getFilter().upper()
            # Keep the same table name
            else:
                tableName = self.tableName

        # Setup the boundary of magnitude based on the filter
        lowMagnitude, highMagnitude = self.filter.getMagBoundary()

        # Get corners of wavefront sensors for this observation field
        obs = self.__getObs(pointing, cameraRotation, self.cameraMJD)

        # Assign the wavefront sensors 
        wavefrontSensors = []
        if (self.camera.name == self.camera.COMCAM):
            if orientation in ("corner", "center", "all"):
                wavefrontSensors = self.camera.getSensor(obs, orientation)
        elif (self.camera.name == self.camera.LSST):
            if (orientation == "corner"):
                wavefrontSensors = self.camera.getWavefrontSensor(obs)
            elif (orientation == "all"):
                wavefrontSensors = self.camera.getScineceSensor(obs)

        if not (wavefrontSensors):
            print("No wavefront sensor is allocated.")

        print("Boresight: (RA, Decl) = (%f, %f) " % (pointing[0], pointing[1]))

        # Query the star database
        starMap = {}
        neighborStarMap = {}
        for detector in wavefrontSensors:

            print("Processing detector %s" % detector)
            wavefrontSensor = wavefrontSensors[detector]

            # Get stars in this wavefront sensor for this observation field
            stars = self.db.query(tableName, cameraFilter, wavefrontSensor[0], wavefrontSensor[1], 
                                    wavefrontSensor[2], wavefrontSensor[3])

            starsQueried = len(stars.RA)
            print("\t\tStars queried: %d" % starsQueried)

            # Populate detector information for the stars
            stars.populateDetector(detector)

            # Populate pixel information for stars
            self.camera.populatePixelFromRADecl(stars, obs)

            # Remove stars that are not on the detector
            self.camera.removeStarsNotOnDetectorSimple(stars, obs, offset)
            starMap[detector] = stars

            starsOnDetector = len(stars.RA)
            print("\t\tStars on detector: %d" % starsOnDetector)

            # Check the candidate of bright stars based on the magnitude
            indexCandidate = stars.checkCandidateStars(cameraFilter, lowMagnitude, highMagnitude)

            # Determine the neighboring stars based on the distance and allowed number 
            # of neighboring stars
            neighborStar = stars.getNeighboringStar(indexCandidate, self.maxDistance, 
                                                    cameraFilter, self.maxNeighboringStar)
            neighborStarMap[detector] = neighborStar

            print("\t\tAvailable candidate stars: %d" % len(neighborStar.SimobjID))

        return neighborStarMap, starMap, wavefrontSensors

    def __getObs(self, pointing, cameraRotation, mjd):
        """
        
        Get the observation metadata.
        
        Arguments:
            pointing {[tuple]} -- Camera boresight (RA, Decl) in degree.
            cameraRotation {[float]} -- Camera rotation angle in degree.
            mjd {[float]} -- Camera mjd.
        
        Returns:
            [metadata] -- The observation meta data (found in the lsst-sims stack) that defines 
                          the pointing.
        """

        # Get Ra, Decl
        RA, Dec = pointing

        obs = []
        obs = ObservationMetaData(pointingRA = RA, pointingDec = Dec, 
                                    rotSkyPos = cameraRotation, mjd = mjd)
        if (not obs):
            print("No observation metadata.")

        return obs

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

    def configNbrCriteria(self, starRadiusInPixel, spacingCoefficient, maxNeighboringStar=99):
        """
        
        Set the neighboring star criteria to decide the scientific target.

        Arguments:
            starRadiusInPixel {[float]} -- Diameter of star. For the defocus = 1.5 mm, the star's radius 
                                           is 63 pixel.
            spacingCoefficient {[float]} -- Maximum distance in units of radius one donut must be 
                                            considered as a neighbor.
        
        Keyword Arguments:
            maxNeighboringStar {[int]} -- Maximum number of neighboring stars. (default: {99})
        """

        self.maxDistance = starRadiusInPixel*spacingCoefficient
        self.maxNeighboringStar = int(maxNeighboringStar)

    def setFilter(self, atype):
        """
        
        Set the active filter type.
        
        Arguments:
            atype {[string]} -- Filter type.
        """

        self.filter.setFilter(atype)

    def getFilter(self):
        """
        
        Get the active filter type.
        
        Returns:
            [string] -- Filter type.
        """

        return self.filter.getFilter()

    def getStddevSplit(self):
        """
        
        Get the standard deviation in database to decide the camera crosses RA=0 or not.
        
        Returns:
            [float] -- Standard deviation.
        """

        return self.db.stddevSplit

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

class SourceSelectorTest(unittest.TestCase):
    """
    Test the source selector. 
    """

    def setUp(self):

        # Camera type: "lsst" or "comcam"
        cameraType = "comcam"

        # Active filter type
        aFilterType = "r"

        # Address of local database
        dbAdress = os.path.join("..", "test", "bsc.db3")

        # Remote database setting
        databaseHost = "localhost:51433"
        databaseUser = "LSST-2"
        databasePassword = "L$$TUser"
        databaseName = "LSSTCATSIM"


        # Set the database
        self.remoteDb = SourceSelector()
        self.localDb = SourceSelector()

        self.remoteDb.configSelector(cameraType=cameraType, dbType="UWdb", aFilter=aFilterType)
        self.localDb.configSelector(cameraType=cameraType, dbType="LocalDb", aFilter=aFilterType)

        # Remote database infomation
        remoteDbInfo = [databaseHost, databaseUser, databasePassword, databaseName]

        # Connect to database
        self.remoteDb.connect(*remoteDbInfo)
        self.localDb.connect(dbAdress)

    def tearDown(self):

        # Disconnect database
        self.remoteDb.disconnect()
        self.localDb.disconnect()

    def testFunctions(self):

        # Boresight (RA, Dec) (unit: degree) (0 <= RA <= 360, -90 <= Dec <= 90)
        pointing = (20.0, 30.0)

        # Camera rotation
        cameraRotation = 0.0

        # Camera orientation for ComCam ("center" or "corner" or "all")
        # Camera orientation for LSSTcam ("corner" or "all")
        orientation = "center"

        # Maximum distance in units of radius one donut must be considered as a neighbor.
        spacingCoefficient = 2.5

        # For the defocus = 1.5 mm, the star's radius is 63 pixel.
        starRadiusInPixel = 63

        # Set the configuration to select the scientific target
        self.remoteDb.configNbrCriteria(starRadiusInPixel, spacingCoefficient, maxNeighboringStar=99)
        self.localDb.configNbrCriteria(starRadiusInPixel, spacingCoefficient, maxNeighboringStar=99)

        # Set the active filter
        # self.remoteDb.setFilter(self.aFilterType)
        # self.localDb.setFilter(self.aFilterType)

        # Test to get the active filter
        self.assertEqual(self.localDb.getFilter(), "r")

        # Test to get the standard deviation 
        self.assertEqual(self.localDb.getStddevSplit(), 20.0)

        # Get the scientific target by querying the remote database
        neighborStarMap, starMap, wavefrontSensors = self.remoteDb.getTargetStar(pointing, 
                                                        cameraRotation, orientation=orientation)

        # Test to get at least one star
        allStars = starMap["R:2,2 S:1,1"]
        self.assertTrue(len(allStars.SimobjID)>=1)

        # Get the scientific target by querying the local database
        neighborStarMapLocal, starMapLocal, wavefrontSensorsLocal = self.localDb.getTargetStar(pointing, 
                                                                    cameraRotation, orientation=orientation)

        # Test the get the empty star map
        allStarsLocal = starMapLocal["R:2,2 S:1,1"]
        self.assertEqual(allStarsLocal.SimobjID, [])

        # Insert the neighboring star map into the database
        self.localDb.insertToBSC(neighborStarMap)

        # Query the local database again
        neighborStarMapLocal, starMapLocal, wavefrontSensorsLocal = self.localDb.getTargetStar(pointing, 
                                                                    cameraRotation, orientation=orientation)

        # Test to get all neighboring stars
        allNeighborStarLocal = neighborStarMapLocal["R:2,2 S:1,1"]
        allNeighborStar = neighborStarMap["R:2,2 S:1,1"]
        self.assertEqual(len(allNeighborStarLocal.SimobjID), len(allNeighborStar.SimobjID))

        # Test to trim the margin
        self.remoteDb.trimMargin(neighborStarMap, 1000)

        # Test to search the id of star based on (ra, decl)
        searchStarId = self.localDb.searchRaDecl(20.088157, 29.983533)

        # Test to update the value
        self.localDb.updateBSC([searchStarId[0][0], searchStarId[0][0]], ["ra", "decl"], [200, 200])
        newSearchStarId = self.localDb.searchRaDecl(200, 200)
        self.assertEqual(searchStarId, newSearchStarId)

        # Delete all data in local database
        allStarList = np.arange(1,len(allNeighborStarLocal.RaDecl)+1)
        self.localDb.db.deleteData(self.localDb.getFilter(), allStarList.tolist())

    def testCoorFun(self):
        fitsFilePath = os.path.join("..", "test", "eimage", "v99999999-fr", "E000", "R22", 
                                    "eimage_99999999_R22_S11_E000.fits.gz")
        raList = [0]
        decList = [0]
        xPosList, yPosList = calcPixPos(fitsFilePath, raList, decList)
        self.assertEqual((xPosList[0], yPosList[0]), (2000, 2036))

if __name__ == "__main__":

    # Do the unit test
    unittest.main()