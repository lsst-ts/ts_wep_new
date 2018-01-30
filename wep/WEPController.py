import os
import numpy as np

from wep.WFDataCollector import WFDataCollector
from wep.EimgIsrWrapper import EimgIsrWrapper
from wep.SourceSelector import SourceSelector
from wep.SourceProcessor import SourceProcessor
from wep.WFEstimator import WFEstimator

from wep.IsrWrapper import getImageData

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

    def importPhoSimDataToButler(self, dataDir, obsId=None, aFilter=None, atype=None, 
                                 overwrite=False):
        """
        
        Import the PhoSim simulated data to match with the data butler to use. This means the 
        registry.sqlite3 repo will be inserted with the meta data if necessary.
        
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

    def getButlerData(self, datasetType, dataId=None, immediate=True):
        """
        
        Retrieves a dataset given an input collection data id.
        
        Arguments:
            datasetType {[str]} -- The type of dataset to retrieve.
        
        Keyword Arguments:
            dataId {[dict]} -- The data id. (default: {None})
            immediate {bool} -- If False use a proxy for delayed loading. (default: {True})
        
        Returns:
            [ExposureU] -- Exposure data.
        """

        return self.dataCollector.butler.get(datasetType=datasetType, dataId=dataId, 
                                             immediate=immediate)

    def getDefocalImg(self, snap, raft, sensor, intraObsId, extraObsId, datasetType="eimage", 
                            immediate=True):
        """
        
        Get the defocal images. The defocal images are based on the visits with different piston 
        position.
        
        Arguments:
            snap {[int]} -- Snap number.
            raft {[str]} -- Raft name (e.g. "2,2").
            sensor {[str]} -- Sensor name (e.g. "1,1").
            intraObsId {[int]} -- Observation ID of intra-focal exposure.
            extraObsId {[int]} -- Observation ID of extra-focal exposure.
        
        Keyword Arguments:
            datasetType {str} -- The type of dataset to retrieve. (default: {"eimage"})
            immediate {bool} -- If False use a proxy for delayed loading. (default: {True})
        
        Returns:
            [ndarray] -- Intra-focal image.
            [ndarray] -- Extra-focal image.
        """

        # Get the images
        intraImg = self.__getImgData(intraObsId, snap, raft, sensor, datasetType=datasetType, 
                                     immediate=immediate)
        extraImg = self.__getImgData(extraObsId, snap, raft, sensor, datasetType=datasetType, 
                                     immediate=immediate)

        return intraImg, extraImg

    def __getImgData(self, obsId, snap, raft, sensor, datasetType, immediate):
        """
        
        Get the image data from exposure by data butler.
        
        Arguments:
            obsId {[int]} -- Observation ID.
            snap {[int]} -- Snap number.
            raft {[str]} -- Raft name (e.g. "2,2").
            sensor {[str]} -- Sensor name (e.g. "1,1").
            datasetType {[str]} -- The type of dataset to retrieve.
            immediate {[bool]} -- If False use a proxy for delayed loading.
        
        Returns:
            [ndarray] -- Image data.
        """
        
        # Enforce the number type
        snap = int(snap)
        obsId = int(obsId)

        # Set the data Id
        dataId = dict(visit=obsId, snap=snap, raft=raft, sensor=sensor)

        # Get the exposure 
        exposure = self.getButlerData(datasetType, dataId=dataId, immediate=immediate)

        # Get the image
        img = getImageData(exposure)

        return img

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
    extraObsId = 9007000
    intraObsId = 9007001
    obsIdList = [extraObsId, intraObsId]
    aFilter = "g"

    dataDirList = ["realComCam/output/Extra", "realComCam/output/Intra"]
    atype = "raw"
    for ii in range(2):
        wepCntlr.importPhoSimDataToButler(dataDirList[ii], obsId=obsIdList[ii], aFilter=aFilter, atype=atype, overwrite=False)

    # Get the image data
    snap = 0
    raft = "2,2"
    sensor = "1,1"
    intraImg, extraImg = wepCntlr.getDefocalImg(snap, raft, sensor, intraObsId, extraObsId, datasetType="eimage")

    # Get the image data

