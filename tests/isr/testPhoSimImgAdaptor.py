import os
import shutil
from astropy.io import fits
import unittest

from lsst.ts.wep.isr.PhoSimImgAdaptor import PhoSimImgAdaptor
from lsst.ts.wep.Utility import getModulePath


class TestPhoSimImgAdaptor(unittest.TestCase):
    """Test the PhoSimImgAdaptor class."""

    def setUp(self):

        # Get the path of module
        modulePath = getModulePath()

        # Path of data folder
        self.dataFolderPath = os.path.join(modulePath, "tests", "testData",
                                           "phosim_out")

        # Destination folder
        self.destinationFolder = os.path.join(modulePath, "tests", "testData",
                                              "phosim_out","isrAnalysis")

    def testFunction(self):

        # Image and calibration products folder
        imageFolder = os.path.join("real", "output")
        biasFolder = os.path.join("biasFrame", "output")
        darkCurrentFolder = os.path.join("darkCurrent", "output")
        flatDomeFolder = os.path.join("flatDome", "output")

        # Declare the ISR image
        phosimImage = PhoSimImgAdaptor(self.dataFolderPath,
                                       self.destinationFolder)
        phosimImage.config(rawDir=imageFolder, biasDir=biasFolder,
                           darkDir=darkCurrentFolder, flatDir=flatDomeFolder)

        # Test to get the header file 
        headerFile = phosimImage.collectHeaderInfo()
        self.assertEqual(list(headerFile), ["a_R22_S11_C01"])

        # Test to rearrange the files (raw and eimage) for butler to use
        obsId = 99999999
        aFilter = "r"
        ampImgName, elecImgName = phosimImage.rearrangeFileForButler(
                    aVisit=obsId, eVisit=obsId, aFilter=aFilter, atype="raw")
        self.assertEqual(ampImgName,
                         ['imsim_99999999_R22_S11_C01_E000.fits.gz'])
        self.assertEqual(elecImgName, []) 

        ampImgNameBias, elecImgNameBias = phosimImage.rearrangeFileForButler(
                                                        aVisit=0, atype="bias")  
        ampImgNameDC, elecImgNameDC = phosimImage.rearrangeFileForButler(
                                                        aVisit=1, atype="dark")
        ampImgNameFlat, elecImgNameFlat = phosimImage.rearrangeFileForButler(
                                    aVisit=2, aFilter=aFilter, atype="flat")

        # Check the existence of file
        filePath = os.path.join(self.destinationFolder, "raw", "v99999999-fr", 
                                "E000", "R22", "S11", ampImgName[0])
        self.assertTrue(os.path.isfile(filePath))

        filePath = os.path.join(self.destinationFolder, "bias", "v0", 
                                "R22", "S11", ampImgNameBias[0])
        self.assertTrue(os.path.isfile(filePath))

        filePath = os.path.join(self.destinationFolder, "dark", "v1", 
                                "R22", "S11", ampImgNameDC[0])
        self.assertTrue(os.path.isfile(filePath))

        filePath = os.path.join(self.destinationFolder, "flat", "v2-fr", 
                                "R22", "S11", ampImgNameFlat[0])
        self.assertTrue(os.path.isfile(filePath))

        # Test to import data into the data butler and do the check of header file 
        # pathDestination = os.path.join(self.dataFolderPath, self.destinationFolder)
        phosimImage.importToButler(phosimImage.pathDestination, ampImgName,
                                   "raw", aFilter=aFilter)
        phosimImage.importToButler(phosimImage.pathDestination, ampImgNameBias,
                                   "bias")
        phosimImage.importToButler(phosimImage.pathDestination, ampImgNameDC,
                                   "dark") 
        phosimImage.importToButler(phosimImage.pathDestination, ampImgNameFlat,
                                   "flat", aFilter=aFilter)

        # Check the updated header file 
        headerFile = fits.getheader(filePath, 0)
        self.assertEqual(headerFile["CTYPE1"], "RA---TAN")
        self.assertEqual(headerFile["RADESYS"], "ICRS")
        self.assertEqual(headerFile["VERSION"], 13357)

        # Remove the generated directory
        shutil.rmtree(phosimImage.pathDestination)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
