import os

from wep.WFDataCollector import WFDataCollector
from wep.EimgIsrWrapper import EimgIsrWrapper
from wep.SourceSelector import SourceSelector
from wep.SourceProcessor import SourceProcessor
from wep.WFEstimator import WFEstimator

class WEPController(object):

    def __init__(self):
        
        self.dataCollector = WFDataCollector()
        self.isrWrapper = None
        self.sourSelc = None
        self.sourProc = None
        self.wfsEsti = None
        self.middleWare = None

    def config(self, dataCollector=None, isrWrapper=None, sourSelc=None, sourProc=None, wfsEsti=None, 
                middleWare=None):
        
        self.__setVar(dataCollector, "dataCollector")
        self.__setVar(isrWrapper, "isrWrapper")
        self.__setVar(sourSelc, "sourSelc")
        self.__setVar(sourProc, "sourProc")
        self.__setVar(wfsEsti, "wfsEsti")
        self.__setVar(middleWare, "middleWare")

    def __setVar(self, value, attrName):
        """
        
        Set the value of attribute.
        
        Arguments:
            value {[obj]} -- New value.
            attrName {[str]} -- Attribute name to set the value.
        """

        if (value is not None):
            setattr(self, attrName, value)

    def configDataCollector(self, pathOfRawData=None, destinationPath=None, dbAdress=None, 
                            butlerInputs=None, butlerOutputs=None):
        """
        
        Do the configuration of data collector.
        
        Keyword Arguments:
            pathOfRawData {[str]} -- Path of raw data. (default: {None})
            destinationPath {[str]} -- Path to the destination. (default: {None})
            dbAdress {[str]} -- Path to the registry.sqlite3 repo.  (default: {None})
            butlerInputs {[str]} -- Butler input directory. (default: {None})
            butlerOutputs {[str]} -- Butlter output directory. (default: {None})
        """

        self.dataCollector.config(pathOfRawData=pathOfRawData, destinationPath=destinationPath, 
                            dbAdress=dbAdress, butlerInputs=butlerInputs, butlerOutputs=butlerOutputs)

    def importPhoSimDataToButler(self, dataDir, obsId=None, aFilter=None, atype=None, overwrite=False):
        """
        
        Import the PhoSim simulated data to match with the data butler to use. This means the registry.sqlite3 
        repo will be inserted with the meta data if necessary.
        
        Arguments:
            dataDir {[str]} -- PhoSim FITS data directory.
        
        Keyword Arguments:
            obsId {[int]} -- Visit/ observation ID. (default: {None})
            aFilter {[str]} -- Filter name (u, g, r, i, z, y). (default: {None})
            atype {[str]} -- Dataset type. (default: {None})
            overwrite {[boolean]} -- Overwrite the existed files or not. (default: {False})
        
        Raises:
            ValueError -- Not allowed type ("raw", "bias", "dark", "flat").
        """

        self.dataCollector.importPhoSimDataToButler(dataDir=dataDir, obsId=obsId, aFilter=aFilter, 
                                                                    atype=atype, overwrite=overwrite)

if __name__ == "__main__":
    
    # Initiate the WEP Controller
    wepCntlr = WEPController()

    # Configure the WFS data collector
    pathOfRawData = "../test/phosimOutput"
    destinationPath = "../test"
    butlerInputs = "../test"
    butlerOutputs = "../test"
    dbAdress = "../test/registry.sqlite3"

    wepCntlr.configDataCollector(pathOfRawData=pathOfRawData, destinationPath=destinationPath, 
                                    dbAdress=dbAdress, butlerInputs=butlerInputs, 
                                    butlerOutputs=butlerOutputs)

    # Import the PhoSim simulated image
    obsId = 9006000
    snap = 0
    raft = "2,2"
    sensor = "1,1"

    aFilter = "g"

    dataDir = "realWfs/output"
    atype = "raw"
    wepCntlr.importPhoSimDataToButler(dataDir, obsId=obsId, aFilter=aFilter, atype=atype, overwrite=False)

    # Let the butler to hold the data
    
