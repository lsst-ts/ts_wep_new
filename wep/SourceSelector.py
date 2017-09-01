from lsst.sims.utils import ObservationMetaData

from bsc.BrightStarDatabase import BrightStarDatabase
from bsc.LocalDatabase import LocalDatabase

from bsc.CameraData import LsstCamera, ComCam 

from bsc.Filter import Filter
from bsc.PlotStarFunc import plotPixel, plotRaDecl

from numpy import nan

import time

class SourceSelector(object):

	def __init__(self, dbType, cameraType):

		if (dbType == "UWdb"):
			self.db = BrightStarDatabase()
		elif (dbType == "LocalDb"):
			self.db = LocalDatabase()
		else:
			raise ValueError("No '%s' database." % dbType)

		self.name = dbType

		if (cameraType == "lsst"):
			self.camera = LsstCamera()
		elif (cameraType == "comcam"):
			self.camera = ComCam()
		else:
			raise ValueError("No '%s' camera." % cameraType)
		self.camera.initializeDetectors()

		self.maxDistance = nan
		self.maxNeighboringStar = nan

		self.filter = Filter()

	def connect(self, *kwargs):
		self.db.connect(*kwargs)

	def disconnect(self):
		self.db.disconnect()

	def getTargetStar(self, tableName, pointing, cameraRotation, cameraMJD, orientation=None, offset=0):

	    # Filter type
	    cameraFilter = self.filter.getFilter()

	    # Get Ra, Decl
	    RA, Dec = pointing

	    # Regenerate the tablen name for local database
	    if (self.name == "LocalDb"):
	    	tableName = tableName + self.filter.getFilter().upper()
	    
	    # Setup the boundary of magnitude based on the filter
	    lowMagnitude, highMagnitude = self.filter.getMagBoundary()

	    # Get corners of wavefront sensors for this observation field
	    obs = []
	    obs = ObservationMetaData(pointingRA = RA, pointingDec = Dec, 
	                              rotSkyPos = cameraRotation, mjd = cameraMJD)
	    if not (obs):
	        print "No observation metadata."

	    # The CCD dimension here is an estimation. Based on LCA-13381, there are three types of sensors.
	    # e2V CCD250: 40.04 mm x 40.96 mm
	    # STA 4400: 20.00 mm x 40.72 mm
	    # STA 3800C: 40.00 mm x 40.72 mm
	    wavefrontSensors = []

	    # Need to update the method to get the senser dimension
	    if (self.camera.name == self.camera.COMCAM):
	        if orientation in ("corner", "center", "all"):
	            wavefrontSensors = self.camera.getSensor(obs, orientation)
	            sensorDimension = [4072, 4000]
	    elif (self.camera.name == self.camera.LSST):
	        wavefrontSensors = self.camera.getWavefrontSensor(obs)
	        sensorDimension = [4072, 4000]

	    if not (wavefrontSensors):
	        print "No wavefront sensor is allocated."

	    print "Boresight: (RA, Decl) = (%f, %f) " % (RA, Dec)

	    # Query the star database
	    starMap = {}
	    neighborStarMap = {}
	    for detector in wavefrontSensors:
	        
	        print "Processing detector %s" % detector
	        wavefrontSensor = wavefrontSensors[detector]

	        # Get stars in this wavefront sensor for this observation field
	        stars = self.db.query(tableName, cameraFilter, wavefrontSensor[0], wavefrontSensor[1], 
	                              wavefrontSensor[2], wavefrontSensor[3])

	        starsQueried = len(stars.RA)
	        print "\t\tStars queried: %d" % starsQueried
	        
	        # Populate detector information for the stars
	        stars.populateDetector(detector)
	         
	        # Populate pixel information for stars
	        self.camera.populatePixelFromRADecl(stars, obs)
	                
	        # Remove stars that are not on the detector
	        self.camera.removeStarsNotOnDetectorSimple(stars, obs, sensorDimension, offset)
	        starMap[detector] = stars
	        
	        starsOnDetector = len(stars.RA)
	        print "\t\tStars on detector: %d" % starsOnDetector
	    
	        # Check the candidate of bright stars based on the magnitude
	        indexCandidate = stars.checkCandidateStars(cameraFilter, lowMagnitude, highMagnitude)

	        # Determine the neighboring stars based on the distance and allowed number 
	        # of neighboring stars
	        neighborStar = stars.getNeighboringStar(indexCandidate, self.maxDistance, 
	                                                cameraFilter, self.maxNeighboringStar)
	        
	        neighborStarMap[detector] = neighborStar

	        print "\t\tAvailable candidate stars: %d" % len(neighborStar.SimobjID)

	    return neighborStarMap, starMap, wavefrontSensors

	def insertToBSC(self):
		# Insert data into local BSC database
		pass

	def generateBSC(self):
		# Generate the bright star catalog.
		pass

	def updateBSC(self):
		# Update the bright star catalog.
		pass

	def analDis(self):
		# Analyze the distribution of target stars.
		pass

	def selWfs(self):
		# Select the wavefront sensors for generating master images.
		pass

	def configSelect(self, maxDistance, maxNeighboringStar=99):
		"""
		
		Set the configuration to decide the scientific target.

		Arguments:
			maxDistance {[float]} -- Maximum distance in units of pixel to decide the neighboring star.
		
		Keyword Arguments:
			maxNeighboringStar {[int]} -- Maximum number of neighboring stars. (default: {99})
		"""

		self.maxDistance = maxDistance
		self.maxNeighboringStar = int(maxNeighboringStar)

	def setFilter(self, atype):
		self.filter.setFilter(atype)

	def getStddevSplit(self):
		return self.db.stddevSplit

if __name__ == "__main__":

    # Remote database setting
    databaseHost = "localhost:51433"
    databaseUser = "LSST-2"
    databasePassword = "L$$TUser"
    databaseName = "LSSTCATSIM"

    # Table name
    tableName = "bright_stars"

    # Local database address
    dbAdress = "/Users/Wolf/bsc.db3"
    localTableName = "BrightStarCatalog"

	# Instantiate the databases
    cameraType = "comcam"
    remoteDb = SourceSelector("UWdb", cameraType)
    localDb = SourceSelector("LocalDb", cameraType)

    # Connect to database
    remoteDbInfo = [databaseHost,databaseUser, databasePassword, databaseName]
    remoteDb.connect(*remoteDbInfo)
    localDb.connect(dbAdress)

    # Boresight (RA, Dec) (unit: degree) (0 <= RA <= 360, -90 <= Dec <= 90)
    pointing = (0.0, 30.0)

    # Camera rotation
    cameraRotation = 0.0
    cameraMJD = 59580.0

    # Maximum distance in units of radius one donut must be considered as a neighbor.
    spacingCoefficient = 2.5

	# For the defocus = 1.5 mm, the star's radius is 63 pixel.
    starRadiusInPixel = 63

	# Set the configuration to select the scientific target
    maxDistance = starRadiusInPixel*spacingCoefficient
    remoteDb.configSelect(maxDistance, maxNeighboringStar=99)
    localDb.configSelect(maxDistance, maxNeighboringStar=99)

	# Set the active filter
    remoteDb.setFilter("u")
    localDb.setFilter("u")

    # Camera orientation for ComCam ("center" or "corner" or "all")
    orientation = "all"

    # Do the query
    t0 = time.time()

    neighborStarMap, starMap, wavefrontSensors = remoteDb.getTargetStar(tableName, pointing, cameraRotation, cameraMJD, orientation=orientation)
    # neighborStarMap, starMap, wavefrontSensors = localDb.getTargetStar(localTableName, pointing, cameraRotation, cameraMJD, orientation=orientation)
    
    t1 = time.time()
    print (t1-t0)

    # Disconnect database
    remoteDb.disconnect()
    localDb.disconnect()

    # Plot the star map in (ra, dec)
    plotRaDecl(wavefrontSensors, starMap, neighborStarMap, localDb.getStddevSplit())





