import os, subprocess, unittest

class SciWFDataCollector(object):

    def __init__(self):
        """
        
        Initialize the SciWFDataCollector class.
        """
        
        self.pathOfRawData = None
        self.destinationPath = None

    def config(self, pathOfRawData=None, destinationPath=None):
        """
        
        Do the configuration.
        
        Keyword Arguments:
            pathOfRawData {[str]} -- Path of raw data. (default: {None})
            destinationPath {[str]} -- Path to the destination. (default: {None})
        """
        
        self.pathOfRawData = pathOfRawData
        self.destinationPath = destinationPath

    def setMapper(self):
        """
        
        Set the camera mapper as LsstSimMapper.
        """

        if (self.destinationPath is not None):

            command = "echo"

            # Set the camera mapper as LsstSimMapper
            argstring = "'lsst.obs.lsstSim.LsstSimMapper' > %s/_mapper" % self.destinationPath
            
            # Execute the command from shell
            runProgram(command, argstring=argstring)

        else:
            print("Can not set _mapper becasue of no destination directory path.")

    def ingestCalibs(self):
        """
        
        Make the faked flat files and ingest them as the calibration files. This step is 
        time-consuming and only needs to do once.
        """
        
        if (self.destinationPath is not None):

            # Make the gain image files
            command = "makeGainImages.py"
            runProgram(command)

            # Ingest the calibration files
            command = "ingestCalibs.py"

            # Set the argument to ingest the calibration files
            argstring = "%s R*.fits --validity 99999 --output %s" % (self.destinationPath, 
                                                                     self.destinationPath)
            runProgram(command, argstring=argstring)

        else:
            print("Can not ingest calibration files becasue of no destination directory path.")

    def ingestSimImages(self, fitsFileArg="lsst_*.fits.gz"):
        """
        
        Import the PhoSim simulated data to match with the data butler to use. This means the 
        registry.sqlite3 repo will be inserted with the meta data if necessary.
        
        Keyword Arguments:
            fitsFileArg {[str]} -- Fits file argument. (default: {"lsst_*.fits.gz"})
        """
        
        if (self.pathOfRawData is not None) and (self.destinationPath is not None):

            # Ingest the simulation images
            command = "ingestSimImages.py"

            # Set the input directory and phosim raw data
            argstring = "%s %s" % (self.destinationPath, os.path.join(self.pathOfRawData, fitsFileArg))

            # Execute the command from shell
            runProgram(command, argstring=argstring)
            
        else:
            print("Can not ingest fits files becasue of no destination directory path.")

def runProgram(command, binDir=None, argstring=None):
    """
    
    Run the program w/o arguments.
    
    Arguments:
        command {[str]} -- Command of application.
    
    Keyword Arguments:
        binDir {[str]} -- Directory of binary application. (default: {None})
        argstring {[str]} -- Arguments of program. (default: {None})
    
    Raises:
        RuntimeError -- There is the error in running the program.
    """

    # Directory of binary application
    if (binDir is not None):
        command = os.path.join(binDir, command)

    # Arguments for the program
    if (argstring is not None):
        command += (" " + argstring)

    # Call the program w/o arguments
    if (subprocess.call(command, shell=True) != 0):
        raise RuntimeError("Error running: %s" % command)

class SciWFDataCollectorTest(unittest.TestCase):

    """ 
    Test the function of SciWFDataCollector.
    """
    
    def setUp(self):
        
        self.tempDir = "../test/temp"

    def testFunction(self):
        
        # Test to run the shell command
        runProgram("mkdir", argstring=self.tempDir)

        # Initiate the collector object
        sciWfDataCollector = SciWFDataCollector()

        # Do the configuration
        sciWfDataCollector.config(destinationPath=self.tempDir)

        # Test to create the _mapper file
        sciWfDataCollector.setMapper()

        # Check the file content
        aFile = open(os.path.join(self.tempDir, "_mapper"), "r")
        content = aFile.read()
        aFile.close()
        answer = "lsst.obs.lsstSim.LsstSimMapper"
        self.assertEqual(content.split()[0], answer)

        # Do the ingestSimImages()

        # Do the ingestCalibs()

        # Remove the file and directory
        argstring = "-rf %s" % self.tempDir
        runProgram("rm", argstring=argstring)

if __name__ == "__main__":

    # Do the unit test
    unittest.main()