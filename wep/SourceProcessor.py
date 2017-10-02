import os
import numpy as np

from deblend.BlendedImageDecorator import BlendedImageDecorator

from SourceSelector import SourceSelector

from matplotlib.colors import LogNorm, SymLogNorm
import matplotlib.pylab as plt

class SourceProcessor(object):

	def __init__(self, sensorName):
		
		self.sensorName = sensorName
		self.sensorDim = None

		self.blendedImageDecorator = BlendedImageDecorator()

	def config(self, sensorName=None, sensorDim=None):

		# Give the sensor name
		if (sensorName is not None):
			self.sensorName = sensorName

		# Give the dimension of sensor
		if (sensorDim is not None):
			self.sensorDim = sensorDim

	def analFieldXY(self):
		# Analyze the field X and Y

		pass

	def migrateDonut(self):
		# Migrate the donut images to reference point for off-axis correction.
		pass

	def generateMasterImg(self):
		# Generate the master images.
		pass

	def removeBg(self):
		# Remove the backgrond noise.
		pass

	def evalSNR(self):
		# Evaluate the SNR of donut.
		pass

	def analDonutImgQual(self):
		# Analyze the donut image quality.
		pass

	def pairDonutImg(self):
		# Pair the intra- and extra-focal donut images.
		pass

	def correctVignette(self):
		# Correct the vignette of donut images (or skip them?).
		pass

	def overlapDonutImg(self):
		# Overlap donut images.
		pass

	def simulateImg(self, imageFolderPath, defocalDis, neighboringStarMapOnSingleSensor, aFilterType):
		"""
		
		Simulate the defocal CCD images with the neighboring star map.
		
		Arguments:
			imageFolderPath {[str]} -- Path to image directory.
			defocalDis {[float]} -- Defocal distance in mm.
			neighboringStarMapOnSingleSensor {[dict]} -- Neighboring star map.
			aFilterType {[string]} -- Active filter type.
		
		Returns:
			[float] -- Simulated intra- and extra-focal images.
		
		Raises:
			ValueError -- No intra-focal image files.
			ValueError -- Numbers of intra- and extra-focal image files are different.
		"""

		# Generate the intra- and extra-focal ccd images
		ccdImgIntra = np.zeros(self.sensorDim)
		ccdImgExtra = np.zeros(self.sensorDim)

		# Redefine the format of defocal distance
		defocalDis = "%.2f" % defocalDis

		# Get all files in the image directory in a sorted order
		fileList = sorted(os.listdir(imageFolderPath))

		# Get the available donut files
		intraFileList = []
		extraFileList = []
		for afile in fileList:

			# Get the file name
			fileName, fileExtension = os.path.splitext(afile)

			# Split the file name for the analysis
			fileNameStr = fileName.split("_")

			# Find the file name with the correct defocal distance
			if (len(fileNameStr) == 3 and fileNameStr[1] == defocalDis):
				
				# Collect the file name based on the defocal type
				if (fileNameStr[-1] == "intra"):
					intraFileList.append(afile)
				elif (fileNameStr[-1] == "extra"):
					extraFileList.append(afile)

		# Get the number of available files
		numFile = len(intraFileList)
		if (numFile == 0):
			raise ValueError("No available donut images.")

		# Check the numbers of intra- and extra-focal images should be the same
		if (numFile != len(extraFileList)):
			raise ValueError("The numbers of intra- and extra-focal images are different.")

		# Get the magnitude of stars
		nameOfMagAttribute = "LSSTMag" + aFilterType.upper()
		starMag = getattr(neighboringStarMapOnSingleSensor, nameOfMagAttribute)

		# Based on the neighboringStarMapOnSingleSensor to reconstruct the image
		for brightStar, neighboringStar in neighboringStarMapOnSingleSensor.SimobjID.items():

			# Generate a random number
			randNum = np.random.randint(0, high=numFile)
			
			# Choose a random donut image from the file
			donutImageIntra = self.__getDonutImgFromFile(imageFolderPath, intraFileList[randNum])
			donutImageExtra = self.__getDonutImgFromFile(imageFolderPath, extraFileList[randNum])

			# Get the bright star magnitude
			magBS = starMag[brightStar]

			# Combine the bright star and neighboring stars. Put the bright star in the first one.
			allStars = neighboringStar[:]
			allStars.insert(0, brightStar)

			# Add the donut image
			for star in allStars:

				# Get the brigtstar pixel x, y
				starX, starY = neighboringStarMapOnSingleSensor.RaDeclInPixel[star]
				magStar = starMag[star]

				# Ratio of magnitude between donuts (If the magnitudes of stars differs by 5, 
				# the brightness differs by 100.)
				# (Magnitude difference shoulbe be >= 1.)
				magDiff = magStar-magBS
				magRatio = 1/100**(magDiff/5.0)

				# Add the donut image
				self.__addDonutImage(magRatio*donutImageIntra, starX, starY, ccdImgIntra)
				self.__addDonutImage(magRatio*donutImageExtra, starX, starY, ccdImgExtra)

		return ccdImgIntra, ccdImgExtra

	def __getDonutImgFromFile(self, imageFolderPath, fileName):
		"""
		
		Read the donut image from the file.
		
		Arguments:
			imageFolderPath {[str]} -- Path to image directory.
			fileName {[str]} -- File name.
		
		Returns:
			[float] -- Image in numpy array.
		"""

		# Get the donut image from the file by the delegation 
		self.blendedImageDecorator.setImg(imageFile=os.path.join(imageFolderPath, fileName))

		return self.blendedImageDecorator.image.copy()


	def __addDonutImage(self, donutImage, starX, starY, ccdImg):
		"""
		
		Add the donut image to simulated CCD image frame.
		
		Arguments:
			donutImage {[float]} -- Image in numpy array.
			starX {[float]} -- Star position in pixel x.
			starY {[float]} -- Star position in pixel y.
			ccdImg {[float]} -- CCD image in numpy array.
		"""

		# Get the dimension of donut image
		d1, d2 = donutImage.shape

		# Get the interger of position to use as the index
		y = int(starY)
		x = int(starX)

		# Add the donut image on the CCD image
		ccdImg[y-int(d1/2):y-int(d1/2)+d1, x-int(d2/2):x-int(d2/2)+d2] += donutImage

def plotCcdImg(ccdImg, name="", scale="log", vmin=None, vmax=None):
	"""
	
	Plot the ccd image.
	
	Arguments:
		ccdImg {[numpy array]} -- Imgae.
	
	Keyword Arguments:
		name {string} -- Image title name (default: {""}).
		scale {string} -- Scale of image map (log or linear) (default: {"log"}).
		vmin {[float]} -- Mininum value to show. This normalizes the luminance data (default: {None}).
		vmax {[float]} -- Maximum value to show. This normalizes the luminance data (default: {None}).		
	"""

	# Change the scale if needed
	if scale not in ("linear", "log"):
		print("No %s scale to choose. Only 'linear' and 'log' scales are allowed." % scale)
		return

	# Decide the norm in imshow for the ploting
	if (scale == "linear"):
		plotNorm = None
	elif (scale == "log"):
		if (ccdImg.min()) < 0:
			plotNorm = SymLogNorm(linthresh=0.03)
		else:
			plotNorm = LogNorm()
	
	# Plot the image
	plt.figure()
	plt.imshow(ccdImg, origin="lower", norm=plotNorm, vmin=vmin, vmax=vmax)
	plt.colorbar()
	plt.title(name)
	plt.show()

if __name__ == '__main__':

	# Define the database and get the neighboring star map
	# Address of local database
	dbAdress = "/Users/Wolf/bsc.db3"

	# Boresight (RA, Dec) (unit: degree) (0 <= RA <= 360, -90 <= Dec <= 90)
	pointing = (20.0, 30.0)

	# Camera rotation
	cameraRotation = 0.0

	# Active filter type
	aFilterType = "u"

	# Camera type: "lsst" or "comcam"
	cameraType = "comcam"

	# Set the camera MJD
	cameraMJD = 59580.0

	# Camera orientation for ComCam ("center" or "corner" or "all")
	# Camera orientation for LSSTcam ("corner" or "all")
	orientation = "center"

	# Maximum distance in units of radius one donut must be considered as a neighbor.
	spacingCoefficient = 2.5

	# For the defocus = 1.5 mm, the star's radius is 63 pixel.
	starRadiusInPixel = 63

	# Get the neighboring star map
	localDb = SourceSelector("LocalDb", cameraType)
	localDb.connect(dbAdress)
	localDb.config(starRadiusInPixel, spacingCoefficient, maxNeighboringStar=99)
	localDb.setFilter(aFilterType)
	neighborStarMapLocal, starMapLocal, wavefrontSensorsLocal = localDb.getTargetStar(pointing, 
																		cameraRotation, orientation=orientation)
	localDb.disconnect()

	# Get the sensor list
	sensorList = wavefrontSensorsLocal.keys()

	# Donut image folder
	imageFolder = "/Users/Wolf/Documents/stash/cwfs_test_images"
	donutImageFolder = "LSST_C_SN26"
	fieldXY = [0, 0]

	# Instantiate a source processor
	sourProc = SourceProcessor(sensorList[0])

	# Give the path to the image folder
	imageFolderPath = os.path.join(imageFolder, donutImageFolder)
	sourProc.config(sensorDim = (4000, 4072))

	# Generate the simulated image
	defocalDis = 1.5
	ccdImgIntra, ccdImgExtra = sourProc.simulateImg(imageFolderPath, defocalDis, neighborStarMapLocal[sensorList[0]], aFilterType)

	# Plot the ccd image
	plotCcdImg(ccdImgIntra, name="Intra focal image", scale="log")
	plotCcdImg(ccdImgExtra, name="Extra focal image", scale="log")

	plotCcdImg(ccdImgIntra, name="Intra focal image", scale="linear")
	plotCcdImg(ccdImgExtra, name="Extra focal image", scale="linear")


 





