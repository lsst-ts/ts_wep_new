import os
import unittest

from lsst.ts.wep.isr.changePhoSimInstrument import updateData, readData
from lsst.ts.wep.Utility import getModulePath


class TestChangePhoSimInstrument(unittest.TestCase):
    """Test the function to update instument segmentation file."""

    def setUp(self):

        # Get the path of module
        modulePath = getModulePath()

        # Path of data folder
        self.folderPath = os.path.join(modulePath, "tests", "testData")
        self.fileName = "segmentation.txt"

    def testReadData(self):

        # Modigy the segmentaion: parallel prescan, serial overscan,
        # serial prescan, parallel overscan (pixel)
        newData = ["1", "2", "3", "4"]
        sameCcdList = ["R44_S10", "R44_S01"]
        updateData(self.folderPath, self.fileName, newData, "readOutDim",
                   sameCcdList=sameCcdList)

        # Read the file and print the new segmentation
        ccdData = readData(self.folderPath, self.fileName, "readOutDim")
        self.assertEqual(ccdData["R44_S10_C17"], ["4", "0", "1", "0"])
        self.assertEqual(ccdData["R43_S10_C17"], newData)

        # Update back to original one
        newData = ["4", "0", "1", "0"]
        updateData(self.folderPath, self.fileName, newData, "readOutDim")


if __name__ == "__main__":

    # Do the unit test
    unittest.main() 
