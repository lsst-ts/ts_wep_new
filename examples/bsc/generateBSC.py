# This script is to generate the local database of bright star catalog from scratch.
# The tables/ schema in local databse should be setup already. 

import numpy as np
import os, time

from lsst.ts.wep.bsc import BrightStarDatabase, CameraData, LocalDatabase, Filter
from getTargetStar import getTargetStar
from lsst.ts.wep.Utility import getModulePath

if __name__ == "__main__":
    
    # Boresight (unit: degree)
    delta = 0.2
    RaArray = np.arange(0,360,delta)
    DecArray = np.append(np.arange(-90, 90, delta) ,90)

    # Camera rotation
    cameraRotation = 0
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
    camera = CameraData.ComCam()
    camera.initializeDetectors()
  
    # Camera orientation (center" or "corner" or "all" or "lsst")
    # For lsst camerea, orientation = "lsst"
    # For ComCam, orientation = "center" or "corner" or "all"
    
    # orientation = "lsst"
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
    
    t0 = time.time()
    # Generate the local database
    # Do the query
    for RA in RaArray:
	    for Dec in DecArray:
	        neighborStarMap, starMap, wavefrontSensors = getTargetStar(brightStarDatabase, tableName, camera, 
	                                                                   RA, Dec, cameraRotation, cameraMJD, 
	                                                                   activeFilter, maxDistance, maxNeighboringStar, 
	                                                                   orientation=orientation, offset=maxDistance)

	        # Write data into the local database
	        for detector in wavefrontSensors:
	            singleNeighborStarMap = neighborStarMap[detector]
	            if (singleNeighborStarMap):
	                bscDatabase.insertData(activeFilter.getFilter(), singleNeighborStarMap)

    t1 = time.time()
    print("Used time is %f sec." % (t1-t0))

    # Disconnect database
    brightStarDatabase.disconnect()
    bscDatabase.disconnect()
