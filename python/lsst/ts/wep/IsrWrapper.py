import os, unittest
import numpy as np

import lsst.daf.persistence as dafPersistence
from lsst.ip.isr import IsrTask
from lsst.ip.isr.assembleCcdTask import AssembleCcdTask

from lsst.ts.wep.SciIsrWrapper import SciIsrWrapper
from lsst.ts.wep.Utility import getModulePath

class IsrWrapper(SciIsrWrapper):

	def __init__(self):
		"""
		
		Initialize the IsrWrapper class.		
		"""
		super(IsrWrapper, self).__init__()
		self.config = None

	def setConfig(self, doBias=True, doBrighterFatter=False, doDark=True, doDefect=True, doFlat=True, 
				  doFringe=True, doLinearize=True, doWrite=False, overscanFitType="MEDIAN"):
		"""
		
		Set up the configuration of image signature removal (ISR).
		
		Keyword Arguments:
			doBias {bool} -- Apply bias frame correction? (default: {True})
			doBrighterFatter {bool} -- Apply the brighter fatter correction? (default: {False})
			doDark {bool} -- Apply dark frame correction? (default: {True})
			doDefect {bool} -- Apply correction for CCD defects, e.g. hot pixels? (default: {True})
			doFlat {bool} -- Apply flat field correction? (default: {True})
			doFringe {bool} -- Apply fringe correction? (default: {True})
			doLinearize {bool} -- Correct for nonlinearity of the detector's response? (default: {True})
			doWrite {bool} -- Persist postISRCCD? (default: {False})
			overscanFitType {string} -- The method for fitting the overscan bias level. (default: {"MEDIAN"})
		"""

		# Configuration of ISR
		config = IsrTask.ConfigClass()

		# Assemble amp-level exposures into a ccd-level exposure? (default: Ture in ConfigClass())
		config.doAssembleCcd = False

		# Assemble amp-level calibration exposures into ccd-level exposure? (default: False in ConfigClass())
		config.doAssembleIsrExposures = False

		# Apply bias frame correction? (default: Ture in ConfigClass())
		config.doBias = doBias

		# Apply the brighter fatter correction? (default: False in ConfigClass())
		config.doBrighterFatter = doBrighterFatter

		# Apply dark frame correction? (default: Ture in ConfigClass())
		config.doDark = doDark

		# Apply correction for CCD defects, e.g. hot pixels? (default: Ture in ConfigClass())
		config.doDefect = doDefect

		# Apply flat field correction? (default: Ture in ConfigClass())
		config.doFlat = doFlat

		# Apply fringe correction? (default: Ture in ConfigClass())
		config.doFringe = doFringe

		# Correct for nonlinearity of the detector's response? (default: Ture in ConfigClass())
		config.doLinearize = doLinearize

		# Persist postISRCCD? (default: Ture in ConfigClass())
		config.doWrite = doWrite

		# Set the type of overscan correction
		config.overscanFitType = overscanFitType

		# Set the configuration
		self.config = config

	def doISR(self, visit, snap, raft, sensor, channel=None, fakeDatasetType=None, 
				outputDatasetType=None):
		"""
		
		Do the instrument signature removal (ISR). (This is just the initial version for real ISR work. 
		For the evaluation reason, the output is single channel. But it should be the ccd in the 
		final.)
		
		Arguments:
			visit {[int]} -- Visit time.
			snap {int} -- Snap time (0 or 1) means first/ second exposure.
			raft {[str]} -- Raft name.
			sensor {[str]} -- Sensor name.

		Keyword Arguments:
			channel {[str]} -- Channel name. (default: {None})
			fakeDatasetType {[str]} -- Use this type of image supported by lsst camera mapper to 
										simulate the post-ISR image. (default: {None})
			outputDatasetType {[str]} -- Output data type supported by lsst camera mapper. 
										(default: {None})
		
		Returns:
			[exposure] -- Exposure image after ISR.
		"""

		# Exposure image after ISR
		postIsrExposure = None

		# Define the data ID and get the raw amplifer image
		# level in dataRef(): The level of dataId at which to reference
		dataId = dict(visit=visit, snap=snap, raft=raft, sensor=sensor)
		if (channel is not None):
			dataId["channel"] = channel

		ampRef = self.butler.dataRef("raw", level=None, dataId=dataId)
		ampExp = ampRef.get("raw")

		# Do the isr task
		lsstIsrTask = IsrTask(config=self.config)
		isrData = lsstIsrTask.readIsrData(ampRef, ampExp)
		postIsrExposure = lsstIsrTask.run(ampExp, **isrData.getDict()).exposure

		# Put the data into the output path
		if (outputDatasetType is not None):
			self.butler.put(postIsrExposure, outputDatasetType, dataId=dataId)

		return postIsrExposure

	def overscanCorrectAllBias(self, visit, snap, raft, sensor, overscanFitType="MEDIAN"):
		"""
		
		Do the overscan correction for bias frame of single sensor. This will replace the bias frame by the
		overscan-corrected one.
		
		Arguments:
			visit {[int]} -- Visit time.
			snap {int} -- Snap time (0 or 1) means first/ second exposure.
			raft {[string]} -- Raft name.
			sensor {[string]} -- Sensor name.
		
		Keyword Arguments:
			overscanFitType {string} -- The method for fitting the overscan bias level. 
										(default: {"MEDIAN"})
		"""

		# Go through 16 amplifiers
		overscanBiasList = {}
		for ii in range(2):
			for jj in range(8):

				# Get the channel name
				channel = str(ii) + "," + str(jj)

				# Do the overscan correction for the single channel
				self.__overscanCorrect(visit, snap, raft, sensor, channel, 
									   overscanFitType=overscanFitType, doWrite=True)

	def __overscanCorrect(self, visit, snap, raft, sensor, channel, overscanFitType="MEDIAN", doWrite=False):
		"""
		
		Do the overscan correction for bias frame of single channel.
		
		Arguments:
			visit {[int]} -- Visit time.
			snap {int} -- Snap time (0 or 1) means first/ second exposure.
			raft {[string]} -- Raft name.
			sensor {[string]} -- Sensor name.
			channel {[string]} -- Channel name.
		
		Keyword Arguments:
			overscanFitType {string} -- The method for fitting the overscan bias level. (default: {"MEDIAN"})
			doWrite {bool} -- Overwrite the file or not. (default: {False})
		
		Returns:
			[Exposure] -- Overscan corrected bias frame.
		"""

		# Degine the data ID and get the bias frame
		dataId = dict(visit=visit, snap=snap, raft=raft, sensor=sensor, channel=channel)
		biasFrame = self.butler.get("bias", dataId=dataId)

		# Define the configuration
		config = self.setConfig(doBias=False, doBrighterFatter=False, doDark=False, doDefect=False, 
								doFlat=False, doFringe=False, doLinearize=False, doWrite=False, 
								overscanFitType=overscanFitType)

		# Correct the overscan
		# Do the isr task
		lsstIsrTask = IsrTask(config=config)
		overscanBias = lsstIsrTask.run(biasFrame).exposure

		# Overwrite into file
		if (doWrite):

			# Put the overscan corrected bias frame
			self.butler.put(overscanBias, "bias", dataId=dataId)

			# Get the file Name
			raftFolder = "R"+"".join(raft.split(","))
			sensorFolder = "S"+"".join(sensor.split(","))
			fileName = "_".join(["imsim", "0", raftFolder, sensorFolder, "C"+"".join(channel.split(","))]) + ".fits.gz"

			# Get the file path
			filePath = os.path.join(self.outputPath, "bias", "v0", raftFolder, sensorFolder, fileName)

			# Check the existence of file
			if (os.path.exists(filePath)):
				os.remove(filePath)

		return overscanBias

	def assembleAmpImg(self, visit, snap, raft, sensor, atype=None, doISR=False, doWrite=False, 
						outputDatasetType="postISRCCD"):
		"""
		
		Assemble amplifier images into single CCD image.
		
		Arguments:
			visit {[int]} -- Number of visit.
			snap {[int]} -- Snap of exposure image. There are two times of exposure in single visit.
							The snap value can be 0 or 1.
			raft {[string]} -- Name of raft.
			sensor {[string]} -- Name of sensor.

		Keyword Arguments:
			atype {[string]} -- Type of data: raw, bias, flat, dark (default: {"None"}).
			doISR {[bool]} -- Do ISR or not? (True or False) (default: {False}).
			doWrite {bool} -- Persist assembled ccd image? (default: {False}) 
			outputDatasetType {[string]} -- Output data type supported by lsst camera mapper. 
						        			(default: {"postISRCCD"})
		
		Returns:
			[butler] -- Data butler of exposure image.
		"""

		# Assemble the amplifier images to get the single CCD image
		# Contruct the dictionary for amplifier images on certain CCD
		ampExp = {}
		# Go through 16 amplifiers
		for ii in range(2):
			for jj in range(8):
				singleChannel = str(ii) + "," + str(jj)

				if (doISR is True and atype == "raw"):
					ampExp[singleChannel] = self.doISR(visit, snap, raft, sensor, channel=singleChannel)
				elif (doISR is False and atype is not None):
					ampExp[singleChannel] = self.getButlerData(visit, snap, raft, sensor, 
															   channel=singleChannel, atype=atype)
				else:
					return

		# Instantiate the AssembleCcdTask
		AssembleCCD = AssembleCcdTask()
		outExposure = AssembleCCD.assembleCcd(ampExp)

		# Persist the assembled ccd image
		if (doWrite):
			
			# Get the data ID
			dataId = dict(visit=visit, snap=snap, raft=raft, sensor=sensor)

			# Put the data into the output path
			self.butler.put(outExposure, outputDatasetType, dataId=dataId)

		return outExposure

class IsrWrapperTest(unittest.TestCase):
	
	"""	
	Test the function of WfsIsrTask.
	"""

	def setUp(self):

		# Get the path of module
		modulePath = getModulePath()

		# Path of data folder
		dataFolderPath = os.path.join(modulePath, "test")
		self.dataFolderPath = dataFolderPath

	def testWfsIsrTask(self):

		# Setting for the ISR correction
		obsId = 99999999
		snap = 0
		raft = "2,2"
		sensor = "1,1"
		channel = "1,4"

		# Instantiate the WFS ISR task
		isrWrapper = IsrWrapper()
		isrWrapper.configWrapper(inputs=self.dataFolderPath, outputs=self.dataFolderPath)
		isrWrapper.configBulter(self.dataFolderPath)

		# Test the function of data butler
		visit = isrWrapper.butler.queryMetadata("raw", ["visit"], dataId={"snap":0})
		self.assertNotEqual(len(visit), 0)

		# Set ISR configuration
		isrWrapper.setConfig(doBias=True, doBrighterFatter=False, doDark=True, 
							 doDefect=True, doFlat=True, doFringe=True, 
							 doLinearize=True, doWrite=False, overscanFitType="MEDIAN")
		
		# Test the setting of configuration
		self.assertEqual(isrWrapper.config.doWrite, False)

		# Test to get the butler data
		butlerDataRaw = isrWrapper.getButlerData(obsId, snap, raft, sensor, channel=channel, atype="raw")
		self.assertEqual(butlerDataRaw.getMetadata().get("GAIN"), 1.83546)

		# Test to assemble the amplifier images to single CCD
		# ccdExp = isrWrapper.assembleAmpImg(obsId, snap, raft, sensor, atype="raw", doISR=False, doWrite=False)
		# self.assertEqual(ccdExp.getDimensions()[0], 4072)
		# self.assertEqual(ccdExp.getDimensions()[1], 4000)

		# Test to do the ISR
		# Before the ISR. Some values are high for the cosmic ray.
		maxRaw = np.max(butlerDataRaw.getMaskedImage().getImage().getArray())
		self.assertGreater(maxRaw, 1000)

		# postIsrExposure = isrWrapper.doISR(obsId, snap, raft, sensor, channel=channel)
		# self.assertEqual(postIsrExposure.getDimensions()[0], 513)
		# self.assertEqual(postIsrExposure.getDimensions()[1], 2001)

		# The cosmic ray shall be corrected.
		# maxIsr = np.max(postIsrExposure.getMaskedImage().getImage().getArray())
		# self.assertLess(maxIsr, 1)

if __name__ == "__main__":

	# Do the unit test
	unittest.main()
