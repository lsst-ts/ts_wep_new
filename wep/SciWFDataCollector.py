import os

from wep.WFDataCollector import WFDataCollector

class SciWFDataCollector(WFDataCollector):

    def __init__(self):
        super(WFDataCollector, self).__init__()

    def setMapperInDir(self):

        if (self.destinationPath is not None):
            command = "echo 'lsst.obs.lsstSim.LsstSimMapper' > %s/_mapper" % self.destinationPath

        print(command)

    def ingestCalib(self):
        pass

    def importPhoSimDataToButler(self):
        pass

if __name__ == "__main__":

    # Define the path
    destinationPath = "/home/ttsai/Document/phosimObsData/input"    

    # Initiate the collector object
    sciWfDataCollector = SciWFDataCollector()

    # Do the configuration
    sciWfDataCollector.config(destinationPath=destinationPath)

    # Set the mapper in the directory
    sciWfDataCollector.setMapperInDir()

