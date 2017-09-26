import time, re
from collections import Iterable

from SALPY_m2ms import SAL_m2ms

class Middleware(object):
	
	def __init__(self, moduleName):
		
		self.__module = __import__("SALPY_"+moduleName)
		self.salMiddleware = getattr(self.__module, "SAL_"+moduleName)()
		self.retTelData = None

		self.moduleName = moduleName

		self.timeOut = -1
		self.timeStart = None

	def shutDownSal(self):
		"""
		
		Shut down the SAL.
		"""

		self.salMiddleware.salShutdown()

	def setTimeOut(self, timeout):
		"""
		
		Set the time out.
		
		Arguments:
			timeOut {number} -- Waiting time for time out.
		"""

		self.timeOut = timeout
		
	def subInfo(self, atype, topic):
		"""
		
		Subscribe the information of specific topic.
		
		Arguments:
			atype {[str]} -- Type of subscription. This can be "event", "command", or "telemetry".
			topic {[str]} -- Topic name.
		
		Raises:
			ValueError -- No type found.
			Warning -- Time out.
			ValueError -- No topic found.
		"""

		# Check the type
		if atype not in ("event", "command", "telemetry"):
			raise ValueError("This is no '%s' type." % atype)

		# Get the current time
		timeNow = time.time()

		# Check the start time. This is only for the first time of  subscription.
		if self.timeStart is None:
			self.timeStart = timeNow

		# Check the time out
		if (timeNow-self.timeStart > self.timeOut):
			raise Warning("Face the time out (%f s) for no new %s data." % (self.timeOut, atype))
		
		# Check the topic exists or not
		if (atype == "telemetry"):
			topic = self.moduleName + "_" + topic
			status = self.salMiddleware.salTelemetrySub(topic)
		elif (atype == "event"):
			topic = self.moduleName + "_logevent_" + topic
			status = self.salMiddleware.salEvent(topic)
		elif (atype == "command"):
			topic = self.moduleName + "_command_" + topic
			status = self.salMiddleware.salProcessor(topic)

		if (status == -1):
			raise ValueError("There is no '%s' topic in %s." % (topic, atype))
			
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
		if (atype == "telemetry"):
			subFuncName = "getNextSample_" + topic.split("_")[-1]
		elif (atype == "event"):
			subFuncName = "getEvent_" + topic.split("_")[-1]
		elif (atype == "command"):
			subFuncName = "acceptCommand_" + topic.split("_")[-1]

		# Get the retrieval data status
		retStatus = self.__getInfo(subFuncName, data)

		# Reset the time if get the new update
		if (retStatus == 0):
			self.timeStart = timeNow

	def __getInfo(self, subFuncName, data):
		"""
		
		Get the information data.
		
		Arguments:
			subFuncName {[str]} -- Function name of subscription.
			data {[dataOfTelemetry]} -- Telemetry data.
		
		Returns:
			[int] -- Status of getting telemetry data.
		"""

		# Get the retrieval data status
		retStatus = getattr(self.salMiddleware, subFuncName)(data)

		# Need to consider the condition of command
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
			print self.retTelData

		return retStatus

	def pubInfo(self, atype, topic, newData):
		"""
		
		Publish the telemetry of specific topic.
		
		Arguments:
			atype {[str]} -- Type of publication. This can be "event", "command", or "telemetry".
			topic {[str]} -- Topic name.
			newData {[dataOfTelemetry]} -- Telemetry data.
		
		Raises:
			ValueError -- No type found.
			ValueError -- No topic found.
		"""

		# Check the type
		if atype not in ("event", "command", "telemetry"):
			raise ValueError("This is no '%s' type." % atype)

		# Check the topic exists or not
		if (atype == "telemetry"):
			topic = self.moduleName + "_" + topic
			status = self.salMiddleware.salTelemetryPub(topic)
		if (atype == "event"):
			topic = self.moduleName + "_logevent_" + topic
			status = self.salMiddleware.salEvent(topic)
		if (atype == "command"):
			topic = self.moduleName + "_command_" + topic
			status = self.salMiddleware.salCommand(topic)

		if (status == -1):
			raise ValueError("There is no '%s' topic in %s." % (topic, atype))

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
		if (atype == "telemetry"):
			pubFuncName = "putSample_" + topic.split("_")[-1]
			getattr(self.salMiddleware, pubFuncName)(data)
		elif (atype == "event"):
			pubFuncName = "logEvent_" + topic.split("_")[-1]
			# Put the "0" here because there is an required interge for input in SAL.
			getattr(self.salMiddleware, pubFuncName)(data, 0)
		elif (atype == "command"):
			pubFuncName = "issueCommand_" + topic.split("_")[-1]
			cmdId = getattr(self.salMiddleware, pubFuncName)(data)

			# Wait for the command to complete, otherwise to abort
			waitFuncName = "waitForCompletion_" + topic.split("_")[-1]
			if (self.timeOut >0):
				getattr(self.salMiddleware, waitFuncName)(cmdId, self.timeOut)
			else:
				# Set the default timeOut as 5 sec
				getattr(self.salMiddleware, waitFuncName)(cmdId, 5)


if __name__ == "__main__":

	# Module name
	moduleName = "m2ms"

	# Declare the SAL middleware
	wepSal = Middleware(moduleName)

	# Set the time out
	timeOut = 15
	wepSal.setTimeOut(timeOut)

	# Set the topic


	# Publish the data for sepecific topic
	# topic name
	topic = "PowerStatus"
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

	# Event topic name
	eventTopic = "SummaryState"

	# Event data
	SummaryStateValue = 2
	priority = 1
	eventData = {"SummaryStateValue": SummaryStateValue,
				 "priority": priority}

	# Command topic name
	commandTopic = "abort"

	# Command data
	state = 2
	commandData = {"state": state}

	# SAL Loop
	sleepTime = 1
	startTime = time.time()
	# for ii in range(3):
	# 	wepSal.pubInfo("telemetry", topic, newData)
	# 	time.sleep(sleepTime)
	# 	print time.time()-startTime
	# 	wepSal.subInfo("telemetry", topic)
	# 	time.sleep(sleepTime)
	# 	print time.time()-startTime

	# for ii in range(10):
	# 	wepSal.pubInfo("event", eventTopic, eventData)
	# 	time.sleep(sleepTime)
	# 	print time.time()-startTime
	# 	wepSal.subInfo("event", eventTopic)
	# 	time.sleep(sleepTime)
	# 	print time.time()-startTime

	for ii in range(3):
		wepSal.pubInfo("command", commandTopic, commandData)
		time.sleep(sleepTime)
		print time.time()-startTime
		# wepSal.subInfo("command", commandTopic)
		# time.sleep(sleepTime)
		# print time.time()-startTime


	# Turn off the SAL
	wepSal.shutDownSal()




