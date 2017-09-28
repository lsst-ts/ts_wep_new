import time

from wep.Middleware import Middleware

from SALPY_tcsWEP import SAL_tcsWEP

if __name__ == "__main__":

	# Set the time out
	timeOut = 15

	# Set the SAL topic
	topic = "WavefrontError"

	# Data in the topic
	sensorID = "R:2,2 S:1,1"
	annularZerikePolynomials = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
	timestamp = time.time()

	newData = {"sensorID": sensorID,
			   "annularZerikePolynomials": annularZerikePolynomials,
			   "timestamp": timestamp}

	# Declare the SAL middleware
	moduleName = "tcsWEP"
	salWepTelPub = Middleware(moduleName)
	salWepTelSub = Middleware(moduleName)

	# Set the time out
	salWepTelSub.setTimeOut(timeOut)

	# Publish the wavefront error
	salWepTelPub.issueTelemetry(topic, newData)

	# Sleep for 1 sec
	time.sleep(1)

	# Subscribe the wavefront error
	salWepTelSub.getTelemetry(topic)





