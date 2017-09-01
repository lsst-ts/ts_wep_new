class WFDataCollector(object):

	def __init__(self):
		pass

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

	def getHexPisPos(self):
		# Get the hexapod piston position to decide intra- and extra-focal images.
		pass

	def importToButler(self):
		# Import the images into data butler.
		pass