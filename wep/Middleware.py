import time, re
from collections import Iterable

from SALPY_m2ms import SAL_m2ms

class Middleware(object):
	
	def __init__(self, moduleName):
		
		self.__module = __import__("SALPY_"+moduleName)
		self.salMiddleware = getattr(self.__module, "SAL_"+moduleName)()
		self.retTelData = None
		
	def subTelemetry(self, topic, updateTime=1, timeOut=-1):
		"""
		
		Subscribe the telemetry of specific topic.
		
		Arguments:
			topic {[str]} -- Topic name.
		
		Keyword Arguments:
			updateTime {number} -- Period time to do the subscription. (default: {1})
			timeOut {number} -- Waiting time for time out. (default: {-1})
		
		Raises:
			ValueError -- No topic found.
		"""
		
		# Check the topic exists or not
		status = self.salMiddleware.salTelemetrySub(topic)
		if (status == -1):
			raise ValueError("There is no '%s' in telemetry subscriber." % topic)
		else:
			print("'%s' subscriber is ready." % topic)

		# Instantiate the data type
		data = getattr(self.__module, topic+"C")()

		# Analyze the attribute names
		telData = {}
		listAttr = dir(data)
		for attribute in listAttr:

			# Use the regular expression to find the public attributes
			# This regex is for the "__"-begin name (private attribute)
			objMatch = re.match(r"^_", attribute)

			# Only take the unmatched attribute
			if (objMatch is None):
				telData[attribute] = None

		self.retTelData = telData

		# Name of subscribe function
		subFuncName = "getNextSample_" + topic.split("_")[-1]

		# Do the time out loop
		if (timeOut>=0):

			# Set the initial condition
			timeStart = time.time()

			# Get the time now
			timeNow = time.time()
			while (timeNow-timeStart <= timeOut):

				# Get the retrieval data status
				retStatus = self.__getTelemetryData(subFuncName, data)

				# Reset the time if get the new update
				if (retStatus == 0):
					timeStart = timeNow
				
				# Assign the update frequency
				time.sleep(updateTime)

				# Get the time now
				timeNow = time.time()

			# Show the time out information
			print("Face the time out (%f s) for no new telemetry data." % timeOut)

		else:

			while (True):

				# Get the retrieval data status
				retStatus = self.__getTelemetryData(subFuncName, data)

				# Assign the update frequency
				time.sleep(updateTime)

		# Turn off the SAL
		# Need to reconsider this part. Maybe need to put the turnOff SAL and while
		# loop in the higher level class.
		wepSal.salMiddleware.salShutdown()

	def __getTelemetryData(self, subFuncName, data):
		"""
		
		Get the telemetry data.
		
		Arguments:
			subFuncName {[str]} -- Function name of subscription.
			data {[dataOfTelemetry]} -- Telemetry data.
		
		Returns:
			[int] -- Status of getting telemetry data.
		"""

		# Get the retrieval data status
		retStatus = getattr(self.salMiddleware, subFuncName)(data)
		if (retStatus == 0):

			# Get the telemetry data
			for aKey in self.retTelData.keys():
				telData = getattr(data, aKey)

				# Check the data is iterable or not
				if (isinstance(telData, Iterable)):
					self.retTelData[aKey] = list(telData)
				else:
					self.retTelData[aKey] = telData

			# Print the information
			print("Get the '%s' telemetry data." % subFuncName)

		return retStatus

	def pubTelemetry(self, topic, newData):
		"""
		
		Publish the telemetry of specific topic.
		
		Arguments:
			topic {[str]} -- Topic name.
			newData {[dataOfTelemetry]} -- Telemetry data.
		
		Raises:
			ValueError -- No topic found.
		"""

		# Check the topic exists or not
		status = self.salMiddleware.salTelemetryPub(topic)
		if (status == -1):
			raise ValueError("There is no '%s' in telemetry publisher." % topic)

		# Instantiate the data type
		data = getattr(self.__module, topic+"C")()

		# Update the data here
		for aKey, aItem in newData.items():

			# Get the needed attribute
			dataItem = getattr(data, aKey)
			
			# Check the item is iterable or not
			if (isinstance((aItem), Iterable)):
				# Put the value one by one
				for ii in range(len(aItem)):
					dataItem[ii] = aItem[ii]
			else:
				# Put the value to data attribute
				setattr(data, aKey, aItem)

		# Publish the data
		pubFuncName = "putSample_" + topic.split("_")[-1]
		getattr(self.salMiddleware, pubFuncName)(data)

if __name__ == "__main__":

	# Module name
	moduleName = "m2ms"

	# Declare the SAL middleware
	wepSal = Middleware(moduleName)

	# Publish the data for sepecific topic
	# topic name
	topic = "m2ms_PowerStatus"
	# topic = "m2ms_MirrorPositionMeasured"

	# Data information of "m2ms_PowerStatus"
	currents = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
	onOff = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
	states = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
	voltages = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]

	newData = {"currents": currents,
			   "onOff": onOff,
			   "states": states,
			   "voltages": voltages}

	# # Data information of "m2ms_MirrorPositionMeasured"
	# xTilt = 0.1
	# yTilt = 0.2
	# piston = 0.3
	# xPosition = 1.0
	# yPosition = 2.0
	# theta_z_position = 3.0

	# newData = {"xTilt": xTilt,
	# 		   "yTilt": yTilt,
	# 		   "piston": piston,
	# 		   "xPosition":xPosition,
	# 		   "yPosition":yPosition,
	# 		   "theta_z_position": theta_z_position}

	# Put the data into a distionary
	wepSal.pubTelemetry(topic, newData)

	# Test to subscribe the data
	updateTime = 1
	timeOut = 15
	# wepSal.subTelemetry(topic, updateTime=updateTime, timeOut=timeOut)
	# wepSal.subTelemetry(topic, updateTime=updateTime)

	# Turn off the SAL
	# wepSal.salMiddleware.salShutdown()




