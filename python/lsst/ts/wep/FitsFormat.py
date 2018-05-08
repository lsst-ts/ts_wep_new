import os, re, time, unittest
import numpy as np
from astropy.io import fits

from lsst.ts.wep.bsc.Filter import Filter
from lsst.ts.wep.Utility import getModulePath

class FitsFormat(object):

    def __init__(self):

        self.fitsDir = None
        self.fitsFilePath = None

    def config(self, fitsDir=None, fitsFilePath=None):
        """
        
        Do the configuration.
        
        Keyword Arguments:
            fitsDir {[str]} -- Directory to save the new FITS file. (default: {None})
            fitsFilePath {[str]} -- FITS fille path. (default: {None})
        """

        self.fitsDir = fitsDir
        self.fitsFilePath = fitsFilePath

    def writeNewFits(self, data, fitsFileName):
        """
        
        Write the new FITS file.
        
        Arguments:
            data {[ndarray]} -- FITS data.
            fitsFileName {[str]} -- FITS file name (e.g. "temp.fits", "temp.fits.gz").
        
        Returns:
            [str] -- New FITS file path. None if the file exists already.
        """

        # Create a PrimaryHDU object to encapsulate the data
        hdu = fits.PrimaryHDU(data)

        # Create a HDUList to contain the newly created primary HDU
        hdul = fits.HDUList([hdu])

        # Write to file
        fitsFilePath = os.path.join(self.fitsDir, fitsFileName)

        try:
            hdul.writeto(fitsFilePath)
        except Exception as OSError:
            fitsFilePath = None
            print(OSError)

        return fitsFilePath

    def updateHeader(self, headerDict):
        """
        
        Update the header in FITS.
        
        Arguments:
            headerDict {[dict]} -- Header in dictionary.
        """

        # Open the fits file
        hdul = fits.open(self.fitsFilePath, mode="update")

        # Get the header
        header = hdul[0].header

        # Add the new date
        for aKey, aItem in headerDict.items():
            header.set(aKey, aItem)

        # Flush the change of header file
        hdul.flush()

        # Close the file
        hdul.close()

    def getMetaDataFromFileName(self, fileName):
        """
        
        Get the metadata (visit, filter, raft, sensor, channel) from the file name.
        
        Arguments:
            fileName {[str]} -- File name.
        
        Returns:
            [dict] -- Metadata.
        
        Raises:
            RuntimeError -- Can not get the meta data.
        """

        # Get the file name
        fileName = os.path.basename(fileName)

        m = re.match(r"\D*_(\d*)_f(\d)_R(\d)(\d)_S(\d)(\d)_C(\d)(\d)", fileName)
        if m is None:
            raise RuntimeError("Cannot get the metadata from file name")

        # Filter dictionary
        filterDict = {"0": Filter.FilterU, 
                      "1": Filter.FilterG, 
                      "2": Filter.FilterR, 
                      "3": Filter.FilterI, 
                      "4": Filter.FilterZ, 
                      "5": Filter.FilterY}

        # Collect the data
        dataDict = {}
        dataDict["OBSID"] = int(m.groups()[0])
        dataDict["FILTER"] = filterDict.get(m.groups()[1])
        dataDict["CCDID"] = "R%s%s_S%s%s" % m.groups()[2:6]
        dataDict["AMPID"] = "C%s%s" % m.groups()[6:8]

        return dataDict

class FitsFormatTest(unittest.TestCase):

    """ 
    Test the function of FitsFormat.
    """

    def setUp(self):

        modulePath = getModulePath()
        self.testDir = os.path.join(modulePath, "test")

    def testFunction(self):

        # Instantiate the fits class
        fitsFormat = FitsFormat()

        # Do the configuration
        fitsFormat.config(fitsDir=self.testDir)
        self.assertEqual(self.testDir, fitsFormat.fitsDir)

        # Define the data
        data = np.random.rand(50,100)*100
        data = data.astype("uint32")

        # Generate the file
        fitsFileName = "temp_1234_f2_R12_S21_C03.fits.gz"
        fitsFilePath = fitsFormat.writeNewFits(data, fitsFileName)
        self.assertTrue(os.path.isfile(fitsFilePath))

        # Configurate the file path
        fitsFormat.config(fitsFilePath=fitsFilePath)

        # Add the metadata based on the file name
        headerDict = fitsFormat.getMetaDataFromFileName(fitsFilePath)
        fitsFormat.updateHeader(headerDict)

        # Get the header
        header = fits.getheader(fitsFilePath)
        self.assertEqual(header["CCDID"], "R12_S21")

        # Remove the file in the final
        os.remove(fitsFilePath)

if __name__ == "__main__":

    # Do the unit test
    unittest.main()