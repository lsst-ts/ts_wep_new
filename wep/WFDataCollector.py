import os, re

from isr.PhoSimImgAdaptor import PhoSimImgAdaptor
from isr.LocalDatabase import LocalDatabase

class WFDataCollector(object):

	def __init__(self, pathOfRawData, destinationPath):
		
		self.pathOfRawData = pathOfRawData
		self.destinationPath = destinationPath
		self.dbAdress = None

	def config(self, pathOfRawData=None, destinationPath=None, dbAdress=None):

		# Set the path to get the raw data
		if (pathOfRawData is not None):
			self.pathOfRawData = pathOfRawData

		# Set the destination path
		if (destinationPath is not None):
			self.destinationPath = destinationPath

		# Set the local database address (registry.sqlite3) for data butler to use
		if (dbAdress is not None):
			self.dbAdress = dbAdress

	def getWfsImg(self):
		# Get the wavefront images.
		pass

	def getCalExp(self):
		# Get the calibration products: bias, flat, dome.
		pass

	def getHistoricWfsImg(self):
		# Get the historic wavefront images.
		pass

	def getWfsMetaData(self):
		# Get the metadata of wavefront images.
		pass

	def importPhoSimDataToButler(self, dataDir, obsId=None, aFilter=None, atype="raw"):

		# Import the images into data butler and check the local database.

		# Check the database has the record or not
		self.__checkRegistryRepo(dataDir)

		return


		# Check the atype
		if atype not in ("raw", "bias", "dark", "flat"):
			raise ValueError("'%s' type is not allowed." % atype)

		# Declare the PhoSim adaptor
		phosimImage = PhoSimImgAdaptor(self.pathOfRawData, self.destinationPath)

		# Configure the path of PhoSimImgAdaptor and rearrange the files
		if (atype == "raw"):
			phosimImage.config(rawDir=dataDir)
			ampImgName, elecImgName = phosimImage.rearrangeFileForButler(aVisit=obsId, eVisit=obsId, aFilter=aFilter, atype=atype) 
		elif (atype == "bias"):
			phosimImage.config(biasDir=dataDir)
			ampImgName, elecImgName = phosimImage.rearrangeFileForButler(aVisit=0, atype=atype) 
		elif (atype == "dark"):
			phosimImage.config(darkDir=dataDir)
			ampImgName, elecImgName = phosimImage.rearrangeFileForButler(aVisit=1, atype=atype) 
		elif (atype == "flat"):
			phosimImage.config(flatDir=dataDir)
			ampImgName, elecImgName = phosimImage.rearrangeFileForButler(aVisit=2, aFilter=aFilter, atype=atype)

		# Import data into the data butler and do the check of header file
		if (atype == "raw"):
			phosimImage.importToButler(self.destinationPath, ampImgName, atype, aFilter=aFilter)
			phosimImage.importToButler(self.destinationPath, elecImgName, "eimage", aFilter=aFilter)
		elif atype in ("bias", "dark"):
			phosimImage.importToButler(self.destinationPath, ampImgName, atype)
		elif (atype == "flat"):
			phosimImage.importToButler(self.destinationPath, ampImgName, atype, aFilter=aFilter)

	def __checkRegistryRepo(self, dataDir):

		global m

		# Analyze the file name in directory
		# Get all files
		fileList = os.listdir(os.path.join(self.pathOfRawData, dataDir))

		# Use the regular expression to analyze the input name
		for fileName in fileList:

			# fileName = "temp_a_99999999_f2_R02_S31_C00_E012.fits.gz"
			
			# Do the file name match
			m = re.match(r"\S*R(\d)(\d)_S(\d)(\d)\S*E(\d*)", fileName)
			
			# Get the information
			raft = "%s,%s" % m.groups()[0:2]
			sensor = "%s,%s" % m.groups()[2:4]
			snap = int(m.groups()[4])

			print(raft, sensor, snap)



		print(fileList)


if __name__ == "__main__":

	# Local database setting
	dbAdress = "../data/registry.sqlite3"

	# Set the path of raw data
	pathOfRawData = "/Users/Wolf/tcs_phosim_output"

	# Destination folder
	destinationPath = "/Users/Wolf/Documents/stash/ts_tcs_wep/data"

	# Instantiate the WFDataCollector
	wfDataCollector = WFDataCollector(pathOfRawData, destinationPath)

	# Import the PhoSim simulated image
	# Rearrange the files (raw and eimage) for butler to use
	dataDir = "real/output"
	atype = "raw"

	obsId = 99999999
	aFilter = "r"

	wfDataCollector.importPhoSimDataToButler(dataDir, obsId=obsId, aFilter=aFilter, atype=atype)


