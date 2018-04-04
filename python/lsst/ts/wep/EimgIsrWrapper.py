import os, shutil, unittest

from lsst.ts.wep.IsrWrapper import IsrWrapper
from lsst.ts.wep.Utility import getModulePath

class EimgIsrWrapper(IsrWrapper):

	def doISR(self, visit, snap, raft, sensor, channel=None, fakeDatasetType="eimage", 
				outputDatasetType="postISRCCD"):
		"""
		
		Do the instrument signature removal. This is just a simulation. The output is the eimage 
		actually.
		
		Arguments:
			visit {[int]} -- Visit time.
			snap {int} -- Snap time (0 or 1) means first/ second exposure.
			raft {[str]} -- Raft name.
			sensor {[str]} -- Sensor name.

		Keyword Arguments:
			channel {[str]} -- Channel name. (default: {None})
			fakeDatasetType {[str]} -- Use this type of image supported by lsst camera mapper 
										to simulate the post-ISR image. (default: {"eimage"})
			outputDatasetType {[str]} -- Output data type supported by lsst camera mapper. 
										(default: {"postISRCCD"})

		Returns:
			[exposure] -- Faked exposure image after ISR.
		"""

		# Exposure image after ISR
		postIsrExposure = None

		# Get the eimage data
		dataId = dict(visit=visit, snap=snap, raft=raft, sensor=sensor)

		# Pretend the eimage is the post ISR CCD
		postIsrExposure = self.butler.get(fakeDatasetType, dataId=dataId)

		# Put the data into the output path
		if (outputDatasetType is not None):
			self.butler.put(postIsrExposure, outputDatasetType, dataId=dataId)

		return postIsrExposure

class EimgIsrWrapperTest(unittest.TestCase):

	"""
    Test the function of EimageWfsIsrTask.
    """

	def setUp(self):

		# Get the path of module
		self.modulePath = getModulePath()

		# Path of data folder
		dataFolderPath = os.path.join(self.modulePath, "test")
		self.dataFolderPath = dataFolderPath

	def testEimageWfsIsrTask(self):

		# Instantiate the WFS ISR task
		eIsrWrapper = EimgIsrWrapper()
		eIsrWrapper.configWrapper(inputs=self.dataFolderPath, outputs=self.dataFolderPath)
		eIsrWrapper.configBulter(self.dataFolderPath, outputs=self.dataFolderPath)

		# Setting for the ISR correction
		obsId = 99999999
		snap = 0
		raft = "2,2"
		sensor = "1,1"

		# Do the ISR
		postIsrExposure = eIsrWrapper.doISR(obsId, snap, raft, sensor, fakeDatasetType="eimage", 
    									   outputDatasetType="postISRCCD")
    	# Check the existence of file
		filePath = os.path.join(self.dataFolderPath, "postISRCCD", "v99999999-fr", "R22", "S11.fits")
		self.assertTrue(os.path.isfile(filePath))

		# Remove the generated directory
		shutil.rmtree(os.path.join(self.dataFolderPath, "postISRCCD"))

if __name__ == "__main__":

	# Do the unit test
	unittest.main()
