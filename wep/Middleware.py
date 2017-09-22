from SALPY_m2ms import SAL_m2ms

class Middleware(object):
	
	def __init__(self):
		
		self.salMiddleware = SAL_m2ms()
		self.__module = __import__("SALPY_m2ms") 
		
	def subTelemetry(self, item):
		pass

	def pubTelemetry(self, item):

		global data

		# Check the item exists or not
		status = self.salMiddleware.salTelemetryPub(item)
		if (status == -1):
			raise ValueError("Theis is no '%s' in telemetry." % item)

		# Instantiate the data type
		data = getattr(self.__module, item+"C")()

		# Update the data here

		# Publish the data
		pubFuncName = "putSample_" + item.split("_")[-1]
		getattr(self.salMiddleware, pubFuncName)(data)

if __name__ == "__main__":

	# Declare the SAL middleware
	wepSal = Middleware()

	# Publish the data for sepecific item
	# item name
	item = "m2ms_PowerStatus"

	# Data information
	currents = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
	onOff = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
	states = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
	volatges = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]

	wepSal.pubTelemetry(item)


