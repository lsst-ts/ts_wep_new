import os, re, unittest, shutil

import lsst.daf.persistence as dafPersistence

from isr.PhoSimImgAdaptor import PhoSimImgAdaptor
from isr.LocalDatabase import LocalDatabase

class WFDataCollector(object):

	def __init__(self):
		"""
		
		Initialize the WFDataCollector class.
		"""
		
		self.pathOfRawData = None
		self.destinationPath = None
		self.db = None
		self.dbAdress = None
		self.butler = None

	def config(self, pathOfRawData=None, destinationPath=None, dbAdress=None, butlerInputs=None, 
				butlerOutputs=None):
		"""
		
		Do the configuration.
		
		Keyword Arguments:
			pathOfRawData {[str]} -- Path of raw data. (default: {None})
			destinationPath {[str]} -- Path to the destination. (default: {None})
			dbAdress {[str]} -- Path to the registry.sqlite3 repo.  (default: {None})
			butlerInputs {[str]} -- Butler input directory. (default: {None})
			butlerOutputs {[str]} -- Butlter output directory. (default: {None})
		"""

		# Set the path to get the raw data
		self.__setVar(pathOfRawData, "pathOfRawData")

		# Set the destination path
		self.__setVar(destinationPath, "destinationPath")

		# Set the local database address (registry.sqlite3) for data butler to use
		if (dbAdress is not None):
			self.db = LocalDatabase()
			self.__setVar(dbAdress, "dbAdress")

		# Instantiate the data butler
		if (butlerInputs is not None) or (butlerOutputs is not None):
			self.butler = dafPersistence.Butler(inputs=butlerInputs, outputs=butlerOutputs)

	def __setVar(self, value, attrName):
		"""
	    
	    Set the value of attribute.
	    
	    Arguments:
	        value {[obj]} -- New value.
	        attrName {[str]} -- Attribute name to set the value.
	    """

		if (value is not None):
			setattr(self, attrName, value)

	def importPhoSimDataToButler(self, dataDir, obsId=None, aFilter=None, atype=None, overwrite=False):
		"""
		
		Import the PhoSim simulated data to match with the data butler to use. This means the registry.sqlite3 
		repo will be inserted with the meta data if necessary.
		
		Arguments:
			dataDir {[str]} -- PhoSim FITS data directory.
		
		Keyword Arguments:
			obsId {[int]} -- Visit/ observation ID. (default: {None})
			aFilter {[str]} -- Filter name (u, g, r, i, z, y). (default: {None})
			atype {[str]} -- Dataset type. (default: {None})
			overwrite {[boolean]} -- Overwrite the existed files or not. (default: {False})
		
		Raises:
			ValueError -- Not allowed type ("raw", "bias", "dark", "flat").
		"""

		# Check the atype
		if atype not in ("raw", "bias", "dark", "flat"):
			raise ValueError("'%s' type is not allowed." % atype)

		# Check the database has the metadata or not. If it is not, the meta data will be put in the repo.
		if (atype == "raw"):
			self.__checkRegistryRepo(dataDir, obsId=obsId, aFilter=aFilter)

		# Declare the PhoSim adaptor
		phosimImage = PhoSimImgAdaptor(self.pathOfRawData, self.destinationPath)

		# Configure the path of PhoSimImgAdaptor and rearrange the files
		if (atype == "raw"):
			phosimImage.config(rawDir=dataDir)
			ampImgName, elecImgName = phosimImage.rearrangeFileForButler(aVisit=obsId, eVisit=obsId, 
													aFilter=aFilter, atype=atype, overwrite=overwrite) 
		elif (atype == "bias"):
			phosimImage.config(biasDir=dataDir)
			ampImgName, elecImgName = phosimImage.rearrangeFileForButler(aVisit=0, atype=atype, 
																		 overwrite=overwrite) 
		elif (atype == "dark"):
			phosimImage.config(darkDir=dataDir)
			ampImgName, elecImgName = phosimImage.rearrangeFileForButler(aVisit=1, atype=atype, 
																		 overwrite=overwrite) 
		elif (atype == "flat"):
			phosimImage.config(flatDir=dataDir)
			ampImgName, elecImgName = phosimImage.rearrangeFileForButler(aVisit=2, aFilter=aFilter, 
																		 atype=atype, overwrite=overwrite)

		# Import data into the data butler and do the check of header file
		if (atype == "raw"):
			phosimImage.importToButler(self.destinationPath, ampImgName, atype, aFilter=aFilter)
			phosimImage.importToButler(self.destinationPath, elecImgName, "eimage", aFilter=aFilter)
		elif atype in ("bias", "dark"):
			phosimImage.importToButler(self.destinationPath, ampImgName, atype)
		elif (atype == "flat"):
			phosimImage.importToButler(self.destinationPath, ampImgName, atype, aFilter=aFilter)

	def __checkRegistryRepo(self, dataDir, obsId=None, aFilter=None):
		"""
		
		Check the survey meta data exists in registry.sqlite3 repo or not. If it is not, the
		meta data will be put in the repo.
		
		Arguments:
			dataDir {[str]} -- FITS data directory.
		
		Keyword Arguments:
			obsId {[int]} -- Visit ID. (default: {None})
			aFilter {[str]} -- Active filter type. (default: {None})
		"""

		# Analyze the file name in directory
		# Get all files
		fileList = os.listdir(os.path.join(self.pathOfRawData, dataDir))

		# Use the regular expression to analyze the input name
		for fileName in fileList:

			# Do the file name match
			m = re.match(r"\S*_R(\d)(\d)_S(\d)(\d)(?:_C(\d)(\d))?\S*E(\d*)", fileName)

			if (m is not None):
			
				# Get the information
				raft = "%s,%s" % m.groups()[0:2]
				sensor = "%s,%s" % m.groups()[2:4]
				channel = None
				if (m.groups()[4] is not None):
					channel = "%s,%s" % m.groups()[4:6]
				snap = int(m.groups()[-1])

				# Set the data Id
				dataId = {"filter": aFilter, "raft": raft, "sensor": sensor, 
						  "snap": snap, "visit": obsId}

				# Check the input image is eimage or amplifier image
				if (channel is not None):
					# Add channel information
					dataId["channel"] = channel
				else:
					# Change the type to eimage
					atype = "eimage"

				# Query the metadata to check the "visit" exists or not. If it returns a empty list,
				# this means the metadata does not exist. 
				queryVisit = self.butler.queryMetadata("raw", "visit", dataId=dataId)

				# Data exists if do not get the empty queryVisit result
				isExist = False
				if (queryVisit):
					isExist = True

				# Put the meta data if not existed.
				if (not isExist):
					# This function should be replaced by the butler API in the final.
					self.__putFakeData(dataId)

	def __putFakeData(self, dataId):
		"""
		
		Put the faked survay meta data into registry.sqlite3 repo. This function should be
		replaced by the butler API in the final.
		
		Arguments:
			dataId {[dict]} -- Data ID.
		"""

		# Connect the local database
		self.db.connect(self.dbAdress)

		# Put the new data
		taiObs = "1994-07-19T06:49:51.543999744"
		skyTile = 30001
		expTime = 15.0

		# insertData(self, visit, snap, raft, sensor, aFilter, taiObs, skyTile, expTime)
		self.db.insertData(dataId["visit"], dataId["snap"], dataId["raft"], dataId["sensor"], 
						   dataId["filter"], taiObs, skyTile, expTime)

		# Disconnect the local database
		self.db.disconnect()

class WFDataCollectorTest(unittest.TestCase):

	"""	
	Test the function of WFDataCollector.
	"""

	def setUp(self):

		# Set the path of raw data
		pathOfRawData = "../test/phosimOutput"

		# Destination folder
		destinationPath = "../test"

		# Butler
		butlerInputs = "../test"
		butlerOutputs = "../test"

		# Instantiate the WFDataCollector
		self.wfDataCollector = WFDataCollector()
		self.wfDataCollector.config(pathOfRawData=pathOfRawData, destinationPath=destinationPath, 
									butlerInputs=butlerInputs, butlerOutputs=butlerOutputs)

	def testFunction(self):

		# Local database setting
		dbAdress = "../test/registry.sqlite3"
		
		# DO the configuration
		self.wfDataCollector.config(dbAdress=dbAdress)
		self.assertEqual(self.wfDataCollector.dbAdress, dbAdress)

		# Import the PhoSim simulated image
		obsId = 99999998
		snap = 0
		raft = "2,2"
		sensor = "1,1"
		aFilter = "r"

		# Search the visit id first
		self.wfDataCollector.db.connect(self.wfDataCollector.dbAdress)
		item = self.wfDataCollector.db.query("raw", visit=obsId, snap=snap, raft=raft, sensor=sensor, aFilter=aFilter)
		self.wfDataCollector.db.disconnect()
		self.assertEqual(item, [])

		# Import the PhoSim simulated image
		dataDir = "real/output"
		atype = "raw"
		self.wfDataCollector.importPhoSimDataToButler(dataDir, obsId=obsId, aFilter=aFilter, atype=atype, overwrite=False)

		# Search the visit id again
		self.wfDataCollector.db.connect(self.wfDataCollector.dbAdress)
		item = self.wfDataCollector.db.query("raw", visit=obsId, snap=snap, raft=raft, sensor=sensor, aFilter=aFilter)
		self.wfDataCollector.db.disconnect()
		self.assertNotEqual(item, [])

		# Delete the meta data in repo
		self.wfDataCollector.db.connect(self.wfDataCollector.dbAdress)
		self.wfDataCollector.db.deleteData(obsId, snap, raft, sensor, aFilter)
		self.wfDataCollector.db.disconnect()

		# Delete the file
		dirName = "v" + str(obsId) + "-f" + str(aFilter)
		dirPath = os.path.join(self.wfDataCollector.destinationPath, atype, dirName)
		shutil.rmtree(dirPath)

if __name__ == "__main__":

	# Do the unit test
	unittest.main()
