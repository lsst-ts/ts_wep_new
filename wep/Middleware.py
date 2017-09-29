import time, re
from collections import Iterable

from SALPY_m2ms import SAL_m2ms

import unittest

class Middleware(object):

	def __init__(self, moduleName):

		self.moduleName = moduleName		
		self.__module = __import__("SALPY_" + moduleName)
		self.salMiddleware = getattr(self.__module, "SAL_" + moduleName)()

		self.retData = None
		self.status = None

		self.timeOut = -1
		self.timeStart = None

	def resetTopic(self):
		"""
		
		Reset the topic.
		"""

		self.retData = None
		self.status = None

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

	def getEvent(self, topic):
		"""
		
		Get the event for specific topic.
		
		Arguments:
			topic {[str]} -- Topic name.
		
		Raises:
			ValueError -- No topic found.
		"""
		
		# Set start time if necessary, and check the timeout
		self.__setStartTimeAndCheckTimeOut()

		# Check the topic exists or not
		salTopic = self.moduleName + "_logevent_" + topic
		if (self.status is None):
			self.status = self.salMiddleware.salEvent(salTopic)

			# Raise the error if no topic is found
			if (self.status == -1):
				raise ValueError("There is no '%s' topic in event." % topic)
			else:
				print("SAL Event: '%s' is ready." % salTopic)

		# Instantiate the data type
		data = getattr(self.__module, salTopic+"C")()

		# Retrieve the data
		if (self.retData is None):
			self.retData = self.__getPubAttr(data)

		# Name of subscribe function
		subFuncName = "getEvent_" + topic

		# Get the retrieval data status
		retStatus = self.__getDataValue(subFuncName, data)

	def getTelemetry(self, topic):
		"""
		
		Get the telemetry for specific topic.
		
		Arguments:
			topic {[str]} -- Topic name.
		
		Raises:
			ValueError -- No topic found.
		"""
		
		# Set start time if necessary, and check the timeout
		self.__setStartTimeAndCheckTimeOut()

		# Check the topic exists or not
		salTopic = self.moduleName + "_" + topic
		if (self.status is None):
			self.status = self.salMiddleware.salTelemetrySub(salTopic)

			# Raise the error if no topic is found
			if (self.status == -1):
				raise ValueError("There is no '%s' topic in telemetry." % topic)
			else:
				print("SAL Telemetry Subscriber: '%s' is ready." % salTopic)

		# Instantiate the data type
		data = getattr(self.__module, salTopic+"C")()

		# Retrieve the data
		if (self.retData is None):
			self.retData = self.__getPubAttr(data)

		# Name of subscribe function
		subFuncName = "getNextSample_" + topic

		# Get the retrieval data status
		retStatus = self.__getDataValue(subFuncName, data)

	def getCommand(self, topic):
		"""
		
		Get the command for specific topic.
		
		Arguments:
			topic {[str]} -- Topic name.
		
		Raises:
			ValueError -- No topic found.
		"""

		# Set start time if necessary, and check the timeout
		self.__setStartTimeAndCheckTimeOut()

		# Check the topic exists or not
		salTopic = self.moduleName + "_command_" + topic
		if (self.status is None):
			self.status = self.salMiddleware.salProcessor(salTopic)

			# Raise the error if no topic is found
			if (self.status == -1):
				raise ValueError("There is no '%s' topic in command." % topic)
			else:
				print("SAL Processor: '%s' is ready." % salTopic)

		# Instantiate the data type
		data = getattr(self.__module, salTopic+"C")()

		# Retrieve the data
		if (self.retData is None):
			self.retData = self.__getPubAttr(data)

		# Name of subscribe function
		subFuncName = "acceptCommand_" + topic

		# Get the retrieval data status
		cmdId = self.__getDataValue(subFuncName, data, showInfo=False)

		# Acknowledge the finish of command
		# Need to add the details how to check the command is done
		isDone = True
		SAL__CMD_COMPLETE = 303
		if (cmdId > 0 and isDone):
			time.sleep(1)
			ackFuncationName = "ackCommand_" + topic
			getattr(self.salMiddleware, ackFuncationName)(cmdId, SAL__CMD_COMPLETE, 0, "Done : OK")

	def issueEvent(self, topic, newData):
		"""
		
		Issue the event for specific topic.
		
		Arguments:
			topic {[str]} -- Topic name.
			newData {[dict]} -- New data for this topic's SAL data instance.
		
		Raises:
			ValueError -- No topic found.
		"""

		# Check the topic exists or not
		salTopic = self.moduleName + "_logevent_" + topic
		if (self.status is None):
			self.status = self.salMiddleware.salEvent(salTopic)

			# Raise the error if no topic is found
			if (self.status == -1):
				raise ValueError("There is no '%s' topic in event." % topic)
			else:
				print("SAL Event: '%s' is ready." % salTopic)

		# Update the data for this topic
		data = self.__setDataValue(salTopic, newData)

		# Publish the data
		pubFuncName = "logEvent_" + topic
		# Put the "0" here because there is an required interge for input in SAL.
		getattr(self.salMiddleware, pubFuncName)(data, 0)

	def issueTelemetry(self, topic, newData):
		"""
		
		Issue the telemetry for specific topic.
		
		Arguments:
			topic {[str]} -- Topic name.
			newData {[dict]} -- New data for this topic's SAL data instance.
		
		Raises:
			ValueError -- No topic found.
		"""

		# Check the topic exists or not
		salTopic = self.moduleName + "_" + topic
		if (self.status is None):
			self.status = self.salMiddleware.salTelemetryPub(salTopic)

			# Raise the error if no topic is found
			if (self.status == -1):
				raise ValueError("There is no '%s' topic in telemetry." % topic)
			else:
				print("SAL Telemetry Publisher: '%s' is ready." % salTopic)

		# Update the data for this topic
		data = self.__setDataValue(salTopic, newData)

		# Publish the data
		pubFuncName = "putSample_" + topic
		getattr(self.salMiddleware, pubFuncName)(data)

	def issueCommand(self, topic, newData, defaultTimeOut=5):
		"""
		
		Issue the command for specific topic.
		
		Arguments:
			topic {[str]} -- Topic name.
			newData {[dict]} -- New data for this topic's SAL data instance.
		
		Keyword Arguments:
			defaultTimeOut {number} -- Default timeout time if it is not set. (default: {5})
		
		Raises:
			ValueError -- No topic found.
		"""

		# Check the topic exists or not
		salTopic = self.moduleName + "_command_" + topic
		if (self.status is None):
			self.status = self.salMiddleware.salCommand(salTopic)

			# Raise the error if no topic is found
			if (self.status == -1):
				raise ValueError("There is no '%s' topic in command." % topic)
			else:
				print("SAL command: '%s' is ready." % salTopic)

		# Update the data for this topic
		data = self.__setDataValue(salTopic, newData)

		# Publish the data
		pubFuncName = "issueCommand_" + topic
		cmdId = getattr(self.salMiddleware, pubFuncName)(data)

		# Wait for the command to complete, otherwise to abort
		waitFuncName = "waitForCompletion_" + topic
		if (self.timeOut >0):
			getattr(self.salMiddleware, waitFuncName)(cmdId, self.timeOut)
		else:
			# Set the default timeOut
			getattr(self.salMiddleware, waitFuncName)(cmdId, defaultTimeOut)

	def __setStartTimeAndCheckTimeOut(self):
		"""
		
		Set the start time and check the time out.
		
		Raises:
			Warning -- Face the time out.
		"""

		# Get the current time
		timeNow = time.time()

		# Check the start time. This is only for the first time of  subscription.
		if self.timeStart is None:
			self.timeStart = timeNow

		# Check the time out
		if (self.timeOut >= 0):
			if (timeNow-self.timeStart > self.timeOut):
				raise Warning("Face the time out (%f s) for no new data." % self.timeOut)

	def __getPubAttr(self, data):
		"""
		
		Get the public attribute in SAL data instance.
		
		Arguments:
			data {[salDataInstance]} -- Instance of SAL data in specific topic.
		
		Returns:
			[dict] -- Attributes of SAL data instance in dictionary.
		"""

		# Analyze the attribute names
		salData = {}
		listAttr = dir(data)
		for attribute in listAttr:

			# Use the regular expression to find the public attributes
			# This regex is for the "__"-begin name (private attribute)
			objMatch = re.match(r"^_", attribute)

			# Only take the unmatched attribute
			if (objMatch is None):
				salData[attribute] = None

		return salData

	def __getDataValue(self, subFuncName, data, showInfo=True):
		"""
		
		Get the information data.
		
		Arguments:
			subFuncName {[str]} -- Function name of subscription.
			data {[dict]} -- New data.
		
		Returns:
			[int] -- Status of getting telemetry data.
		"""

		# Get the retrieval data status
		retStatus = getattr(self.salMiddleware, subFuncName)(data)

		# Need to consider the condition of command
		if (retStatus >= 0):

			# Get the telemetry data
			for aKey in self.retData.keys():
				telData = getattr(data, aKey)

				# Check the data is iterable or not
				# Because the string is iterable in Python, need to consider it also.
				if (isinstance(telData, Iterable) and not isinstance(telData, str)):
					self.retData[aKey] = list(telData)
				else:
					self.retData[aKey] = telData

			# Print the information
			if (showInfo):
				print("Get the '%s' data." % subFuncName)

			# Reset the time if get the new update
			self.timeOut = time.time()

		# Return the status to retrieve the SAL data
		return retStatus

	def __setDataValue(self, salTopic, newData):
		"""
		
		Set the data values in specific sal topic.
		
		Arguments:
			salTopic {[str]} -- SAL topic name.
			newData {[dict]} -- New data in dictionary.
		
		Returns:
			[salDataInstance] -- Instance of SAL data in specific topic.
		"""

		# Instantiate the data type
		data = getattr(self.__module, salTopic+"C")()

		# Update the data here
		for aKey, aItem in newData.items():

			# Get the needed attribute
			dataItem = getattr(data, aKey)
			
			# Check the item is iterable or not
			# Because the string is iterable in Python, need to consider it also.
			if (isinstance(aItem, Iterable) and not isinstance(aItem, str)):
				# Put the value one by one
				for ii in range(len(aItem)):
					dataItem[ii] = aItem[ii]
			else:
				# Put the value to data attribute
				setattr(data, aKey, aItem)

		# Return the sal data instance
		return data

class MiddlewareTest(unittest.TestCase):
	"""
	Test functions in Middleware.
	"""

	def setUp(self):

		# Module name
		moduleName = "m2ms"

		# Declare the Middleware
		self.wepSalIssue = Middleware(moduleName)
		self.wepSalGet = Middleware(moduleName)

	def testTelemetry(self):

		# Reset the topic
		self.wepSalIssue.resetTopic()
		self.wepSalGet.resetTopic()

		# Set the time out
		timeOut = 15
		self.wepSalGet.setTimeOut(timeOut)
		
		# Test the set timeOut
		self.assertEqual(self.wepSalGet.timeOut, timeOut)

		# Set the telemetry topic
		topic = "PowerStatus"

		# Data information of "m2ms_PowerStatus"
		currents = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
		onOff = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
		states = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
		voltages = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]

		newData = {"currents": currents,
				   "onOff": onOff,
				   "states": states,
				   "voltages": voltages}

		# Issue the telemetry
		self.wepSalIssue.issueTelemetry(topic, newData)

		# Sleep 1 sec
		time.sleep(1)

		# Get the telemetry
		self.wepSalGet.getTelemetry(topic)

		# Test to get the telemetry
		self.assertEqual(self.wepSalGet.retData["currents"], currents)

	def testEvent(self):

		# Reset the topic
		self.wepSalIssue.resetTopic()
		self.wepSalGet.resetTopic()

		# Set the event topic
		topic = "SummaryState"

		# Event data
		SummaryStateValue = 2
		priority = 1
		newData = {"SummaryStateValue": SummaryStateValue,
				   "priority": priority}

		# Issue the event
		self.wepSalIssue.issueEvent(topic, newData)

		# Sleep 1 sec
		time.sleep(1)

		# Get the event
		self.wepSalGet.getEvent(topic)

		# Test to get the event information
		self.assertEqual(self.wepSalGet.retData["SummaryStateValue"], SummaryStateValue)

	def testCommand(self):

		# Reset the topic
		self.wepSalIssue.resetTopic()
		self.wepSalGet.resetTopic()

		# Set the command topic
		topic = "abort"

		# Command data
		state = 2
		newData = {"state": state}

		# Issue the event
		self.wepSalIssue.issueCommand(topic, newData)

		# Sleep 1 sec
		self.wepSalGet.getCommand(topic)

		# Test to get the command item list
		self.assertEqual(len(self.wepSalGet.retData), 5)

	def tearDown(self):

		# Turn off the sal
		self.wepSalIssue.shutDownSal()
		self.wepSalGet.shutDownSal()

if __name__ == "__main__":

	# Do the unit test
	unittest.main()
