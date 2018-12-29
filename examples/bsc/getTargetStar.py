# -*- coding: utf-8 -*-
import os, time

from lsst.sims.utils import ObservationMetaData
from lsst.ts.wep.bsc import BrightStarDatabase, CameraData, LocalDatabase, Filter
from lsst.ts.wep.bsc.PlotStarFunc import plotPixel, plotRaDecl
from lsst.ts.wep.Utility import getModulePath

def getTargetStar(database, tableName, camera, RA, Dec, cameraRotation, cameraMJD, 
                  activeFilter, maxDistance, maxNeighboringStar, orientation=None, offset=0):
    """
    
    Analyze the target stars used in source selection of wavefront estimation pipeline (WEP).
    The selection criteria is the magnitude and maximum distance. 
    
    Arguments:
        database {[metadata]} -- Star database.
        tableName {[string]} -- Table name in database.
        camera {[metadata]} -- Camera type. The geometry of comcam is based on the central raft of lsst. 
        RA {[float]} -- RA of telescope boresight.
        Dec {[float]} -- Dec of telescope boresight.
        cameraRotation {[float]} -- Camera rotation angle.
        cameraMJD {[float]} -- Camera MJD.
        activeFilter {[Filter]} -- Active filter.
        maxDistance {[float]} -- Maximum distance in pixel.
        maxNeighboringStar {[float]} -- Maximum number of neighboring stars.

    Keyword Arguments: 
        orientation {[string]} -- Orientation of wavefront sensor(s) on camera. (default: {None})
        offset{[float]} -- Add offset in pixel to sensor dimension for judging stars on detector or not. 
                           offset=0 for normal use. 
                           offset=maxDistance to generate the local database. (default: {0})
        
    Returns:
        neighborStarMap{[list]} -- Information of neighboring stars and candidate stars with 
                                   the name of sensor as a list.
        starMap{[list]} -- Information of stars with the name of sensor as a list.
        wavefrontSensors{[list]} -- Corners of sensor with the name of sensor as a list.
    """

    # Filter type
    cameraFilter = activeFilter.getFilter()
    
    # Setup the boundary of magnitude based on the filter
    lowMagnitude, highMagnitude = activeFilter.getMagBoundary()

    # Get corners of wavefront sensors for this observation field
    obs = []
    obs = ObservationMetaData(pointingRA = RA, pointingDec = Dec, 
                              rotSkyPos = cameraRotation, mjd = cameraMJD)
    if not (obs):
        print("No observation metadata.")

    # Assign the wavefront sensors 
    wavefrontSensors = []
    if (camera.name == camera.COMCAM):
        if orientation in ("corner", "center", "all"):
            wavefrontSensors = camera.getSensor(obs, orientation)
    elif (camera.name == camera.LSST):
        wavefrontSensors = camera.getWavefrontSensor(obs)

    if not (wavefrontSensors):
        print("No wavefront sensor is allocated.")

    print("Boresight: (RA, Decl) = (%f, %f) " % (RA, Dec))

    # Query the star database
    starMap = {}
    neighborStarMap = {}
    for detector in wavefrontSensors:
        
        print("Processing detector %s" % detector)
        wavefrontSensor = wavefrontSensors[detector]

        # Get stars in this wavefront sensor for this observation field
        stars = database.query(tableName, cameraFilter, wavefrontSensor[0], wavefrontSensor[1], 
                               wavefrontSensor[2], wavefrontSensor[3])

        starsQueried = len(stars.RA)
        print("\t\tStars queried: %d" % starsQueried)
        
        # Set the detector information for the stars
        stars.setDetector(detector)
         
        # Populate pixel information for stars
        camera.populatePixelFromRADecl(stars, obs)
                
        # Remove stars that are not on the detector
        camera.removeStarsNotOnDetectorSimple(stars, offset)
        starMap[detector] = stars
        
        starsOnDetector = len(stars.RA)
        print("\t\tStars on detector: %d" % starsOnDetector)
    
        # Check the candidate of bright stars based on the magnitude
        indexCandidate = stars.checkCandidateStars(cameraFilter, lowMagnitude, highMagnitude)

        # Determine the neighboring stars based on the distance and allowed number 
        # of neighboring stars
        neighborStar = stars.getNeighboringStar(indexCandidate, maxDistance, 
                                                cameraFilter, maxNeighboringStar)
        
        neighborStarMap[detector] = neighborStar

        print("\t\tAvailable candidate stars: %d" % len(neighborStar.SimobjID))

    return neighborStarMap, starMap, wavefrontSensors
    
if __name__ == "__main__":
    
    # Boresight (unit: degree)
    RA = 20.0    # 0 <= RA <= 360
    Dec = 30.0  # -90 <= Dec <= 90

    # Camera rotation
    cameraRotation = 0.0
    cameraMJD = 59580.0

    # Maximum distance in units of radius one donut must be to another 
    # to be considered as a neighbor
    spacingCoefficient = 2.5     
    
    # Maximum number of neighboring stars 
    maxNeighboringStar = 99
    
    # Camera information and criteria of star magnitude
    # For the defocus = 1.5 mm, the star's radius is 63 pixel.
    starRadiusInPixel = 63
    activeFilter = Filter.Filter()
    activeFilter.setFilter("u")

    maxDistance = starRadiusInPixel*spacingCoefficient

    # Camera instance (Two classes: LsstCamera and ComCam) 
    # camera = CameraData.LsstCamera()
    camera = CameraData.ComCam()
    camera.initializeDetectors()
  
    # Camera orientation ("center" or "corner" or "all")
    # For ComCam, orientation = "center" or "corner" or "all"
    orientation = "center"
   
    # Remote database setting
    databaseHost = "localhost:51433"
    databaseUser = "LSST-2"
    databasePassword = "L$$TUser"
    databaseName = "LSSTCATSIM"

    # Table name
    tableName = "bright_stars"

    # Get the path of module
    modulePath = getModulePath()

    # Local database address
    dbAdress = os.path.join(modulePath, "test", "bsc.db3")
    localTableName = "BrightStarCatalog" + activeFilter.getFilter().upper()

    # Setup UW database
    brightStarDatabase = BrightStarDatabase.BrightStarDatabase()

    # Setup Local database
    bscDatabase = LocalDatabase.LocalDatabase()

    # Connect to database
    brightStarDatabase.connect(databaseHost,databaseUser, databasePassword, databaseName)
    bscDatabase.connect(dbAdress)
    
    # Do the query
    t0 = time.time()
    neighborStarMap, starMap, wavefrontSensors = getTargetStar(brightStarDatabase, tableName, camera, 
                                                               RA, Dec, cameraRotation, cameraMJD, 
                                                               activeFilter, maxDistance, maxNeighboringStar, 
                                                               orientation=orientation)
    t1 = time.time()
    print(t1-t0)
    # Write data into the local database
    for detector in wavefrontSensors:
        singleNeighborStarMap = neighborStarMap[detector]
        if (singleNeighborStarMap):
            bscDatabase.insertData(activeFilter.getFilter(), singleNeighborStarMap)

    # Disconnect database
    brightStarDatabase.disconnect()
    bscDatabase.disconnect()

    # # Local Database
    # # Setup database
    # bscDatabase = LocalDatabase.LocalDatabase()

    # # Connect to database
    # bscDatabase.connect(dbAdress)

    # # # Delete all data
    # # listID = bscDatabase.getAllId(cameraFilter)
    # # bscDatabase.deleteData(cameraFilter, listID)

    # # Do the query
    # neighborStarMap, starMap, wavefrontSensors = getTargetStar(bscDatabase, localTableName, camera, 
    #                                                            RA, Dec, cameraRotation, cameraMJD, 
    #                                                            cameraFilter, maxDistance, maxNeighboringStar, 
    #                                                            orientation=orientation)
    # # Disconnect database
    # bscDatabase.disconnect()

    # Plot the star map in (ra, dec)
    # plotRaDecl(wavefrontSensors, starMap, neighborStarMap, brightStarDatabase.stddevSplit)

    # Plot the star map in pixel
    # sensor = "R:2,2 S:1,1"
    # plotPixel(starMap[sensor], neighborStarMap[sensor])
    # singleNeighborStarMap = neighborStarMap[sensor]
   
    


 



