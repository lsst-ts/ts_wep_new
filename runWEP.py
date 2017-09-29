import time
import os

from wep.Middleware import Middleware
from wep.WFEstimator import WFEstimator

from SALPY_tcsWEP import SAL_tcsWEP

if __name__ == "__main__":

	# Define the instrument folder
	instruFolder = "/home/te-wei/Documents/stash/ts_lsst_wep_27/instruData"

	# Define the algorithm folder
	algoFolderPath = "/home/te-wei/Documents/stash/ts_lsst_wep_27/algo"

	# Define the image folder and image names
	# Image data -- Don't know the final image format.
	# It is noted that image.readFile inuts is based on the txt file.
	imageFolderPath = "/home/te-wei/Documents/stash/ts_lsst_wep_27/tests/testImages/LSST_NE_SN25"
	intra_image_name = "z11_0.25_intra.txt"
	extra_image_name = "z11_0.25_extra.txt"

	# Path to image files
	intraImgFile = os.path.join(imageFolderPath, intra_image_name)
	extraImgFile = os.path.join(imageFolderPath, extra_image_name)

	# Field XY position
	fieldXY = [1.185, 1.185]

	# Set the time out
	timeOut = 15

	# Set the module name
	moduleName = "tcsWEP"

	# Set the SAL topic
	topic = "WavefrontError"

	# Decalre the WFEsitmator
	wfsEst = WFEstimator(instruFolder, algoFolderPath)

	# Setup the images
	wfsEst.setImg(fieldXY, imageFile=intraImgFile, defocalType="intra")
	wfsEst.setImg(fieldXY, imageFile=extraImgFile, defocalType="extra")

	# Setup the configuration
	# If the configuration is reset, the images are needed to be set again.
	wfsEst.config(solver="exp", debugLevel=0)

	# Declare the SAL middleware
	salWepTelPub = Middleware(moduleName)
	salWepTelSub = Middleware(moduleName)

	# Set the time out
	salWepTelSub.setTimeOut(timeOut)

	# Data in the topic
	sensorID = "R:2,2 S:1,1"
	annularZerikePolynomials = wfsEst.calWfsErr().tolist()
	timestamp = time.time()

	newData = {"sensorID": sensorID,
			   "annularZerikePolynomials": annularZerikePolynomials,
			   "timestamp": timestamp}

	# Publish the wavefront error
	salWepTelPub.issueTelemetry(topic, newData)

	# Sleep for 1 sec
	time.sleep(1)

	# Subscribe the wavefront error
	salWepTelSub.getTelemetry(topic)

	# Turn off the SAL
	salWepTelPub.shutDownSal()
	salWepTelSub.shutDownSal()

	# Print the wavefront error
	print(salWepTelSub.retData["annularZerikePolynomials"])





