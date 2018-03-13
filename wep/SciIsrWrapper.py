import os, unittest

from wep.SciWFDataCollector import runProgram

class SciIsrWrapper(object):

    def __init__(self):
        """
        
        Initialize the SciIsrWrapper class.        
        """

        self.pathData = None
        self.outputPath = None

    def configWrapper(self, inputs=None, outputs=None):
        """
        
        Configurate the ISR wrapper.

        Keyword Arguments:
            inputs {[str]} -- Input directory. (default: {None})
            outputs {[str]} -- Output directory. (default: {None})
        """
        
        self.pathData = inputs
        self.outputPath = outputs

    def doISR(self, visit, snap, raft, sensor):
        """
        
        Do the instrument signature removal (ISR). 
        
        Arguments:
            visit {[int]} -- Visit time.
            snap {int} -- Snap time (0 or 1) means first/ second exposure.
            raft {[str]} -- Raft name.
            sensor {[str]} -- Sensor name.
        """
        
        if (self.pathData is not None) and (self.outputPath is not None):

            # Process the phoSim CCD images
            command = "processSimCcd.py"

            # Set data Id
            dataId = "visit=%d raft=%s sensor=%s snap=%d" % (visit, raft, sensor, snap)

            # Set the input and output directories
            argstring = "%s --id %s --output %s" % (self.pathData, dataId, self.outputPath)

            # Do the ISR and assemble CCD image
            runProgram(command, argstring=argstring)

class SciIsrWrapperTest(unittest.TestCase):


    """ 
    Test the function of SciIsrWrapper.
    """

    def setUp(self):

        # Path of data folder
        dataFolderPath = "../test"
        self.pathData = dataFolderPath
        self.outputPath = dataFolderPath

    def testFunction(self):
        
        # Instantiate the isr wrapper
        isrWrapper = SciIsrWrapper()

        # Configurate the path
        isrWrapper.configWrapper(inputs=self.pathData, outputs=self.outputPath)

        self.assertEqual(isrWrapper.pathData, self.pathData)
        self.assertEqual(isrWrapper.outputPath, self.outputPath)

if __name__ == "__main__":

    # Do the unit test
    unittest.main()
    
    # # Path
    # pathData = "/home/ttsai/Document/phosimObsData/input"
    # outputPath = "/home/ttsai/Document/phosimObsData/output"

    # # Instantiate the isr wrapper
    # isrWrapper = SciIsrWrapper()

    # # Configurate the path
    # isrWrapper.configWrapper(inputs=pathData, outputs=outputPath)

    # # Do the isr
    # visit = 9006000
    # snap = 0
    # raft = "2,2"
    # sensor = "1,1"
    # isrWrapper.doISR(visit, snap, raft, sensor)