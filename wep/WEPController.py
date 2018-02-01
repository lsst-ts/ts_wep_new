import os, re
import numpy as np

from astropy.io import fits

from wep.WFDataCollector import WFDataCollector
from wep.IsrWrapper import getImageData
from wep.EimgIsrWrapper import EimgIsrWrapper
from wep.SourceSelector import SourceSelector
from wep.SourceProcessor import SourceProcessor
from wep.WFEstimator import WFEstimator

from wep.LocalDatabaseDecorator import LocalDatabaseDecorator

class WEPController(object):

    def __init__(self):
        
        self.sourSelc = None
        self.dataCollector = None
        self.isrWrapper = None
        self.sourProc = None
        self.wfsEsti = None
        self.middleWare = None

    def config(self, sourProc=None, dataCollector=None, isrWrapper=None, sourSelc=None, wfsEsti=None, 
                middleWare=None):

        self.__setVar(sourSelc, "sourSelc")        
        self.__setVar(dataCollector, "dataCollector")
        self.__setVar(isrWrapper, "isrWrapper")
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

    def getTargetStarByFile(self, dbAdress, skyInfoFilePath, pointing, cameraRotation, orientation=None, 
                            tableName="TempTable"):
        """
        
        Get the target stars by querying the file.
        
        Arguments:
            dbAdress {[str]} -- Local database address.
            skyInfoFilePath {[str]} -- File path of sky information.
            pointing {[tuple]} -- Camera boresight (RA, Decl) in degree.
            cameraRotation {[float]} -- Camera rotation angle in degree.
        
        Keyword Arguments:
            orientation {[str]} -- Orientation of wavefront sensor(s) on camera. (default: {None})
            tableName {[str]} -- Table name. (default: {None})
        
        Returns:
            neighborStarMap {[list]} -- Information of neighboring stars and candidate stars with 
                                        the name of sensor as a list.
            starMap {[list]} -- Information of stars with the name of sensor as a list.
            wavefrontSensors {[list]} -- Corners of sensor with the name of sensor as a list.
        """

        # Check the database name is local database
        if (self.sourSelc.name != self.sourSelc.LocalDb):
            raise TypeError("The database type is not LocalDatabaseDecorator.")

        # Get the filter type
        aFilter = self.sourSelc.getFilter()

        # Connect the database
        self.sourSelc.connect(dbAdress)

        # Create the table
        self.sourSelc.db.createTable(aFilter, tableName)

        # Insert the sky data
        self.sourSelc.db.insertDataByFile(aFilter, tableName, skyInfoFilePath, skiprows=1)
        
        # Do the query and analysis
        neighborStarMap, starMap, wavefrontSensors = self.sourSelc.getTargetStar(pointing, cameraRotation, 
                                                                orientation=orientation, tableName=tableName)

        neighborStarMap, starMap, wavefrontSensors = self.__analyzeStarMap(neighborStarMap, starMap, 
                                                                                    wavefrontSensors)

        # Delete the table
        self.sourSelc.db.deleteTable(tableName)
        
        # Disconnect the database
        self.sourSelc.disconnect()

        return neighborStarMap, starMap, wavefrontSensors

    def __analyzeStarMap(self, neighborStarMap, starMap, wavefrontSensors):
        """
        
        Analyze the star map and remove the sensor without bright stars.
        
        Arguments:
            neighborStarMap {[list]} -- Information of neighboring stars and candidate stars with 
                                        the name of sensor as a list.
            starMap {[list]} -- Information of stars with the name of sensor as a list.
            wavefrontSensors {[list]} -- Corners of sensor with the name of sensor as a list.
        
        Returns:
            neighborStarMap {[list]} -- Information of neighboring stars and candidate stars with 
                                        the name of sensor as a list.
            starMap {[list]} -- Information of stars with the name of sensor as a list.
            wavefrontSensors {[list]} -- Corners of sensor with the name of sensor as a list.
        """

        # Collect the sensor list without the bright star
        noStarSensorList = []
        for aKey, aItem in starMap.items():
            if len(aItem.RA) == 0:
                noStarSensorList.append(aKey)

        # Remove the keys in map
        for aKey in noStarSensorList:
            neighborStarMap.pop(aKey)
            starMap.pop(aKey)
            wavefrontSensors.pop(aKey)
        
        return neighborStarMap, starMap, wavefrontSensors

    def importPhoSimDataToButler(self, dataDir, atype="raw", overwrite=False):
        """
        
        Import the PhoSim simulated data to match with the data butler to use. This means the 
        registry.sqlite3 repo will be inserted with the meta data if necessary.
        
        Arguments:
            dataDir {[str]} -- PhoSim FITS data directory.
        
        Keyword Arguments:
            atype {[str]} -- Dataset type. (default: {"raw"})
            overwrite {[boolean]} -- Overwrite the existed files or not. (default: {False})
        
        Raises:
            ValueError -- Not allowed type ("raw", "bias", "dark", "flat").
        """

        # Get all files in the directory
        fullDataDir = os.path.join(self.dataCollector.pathOfRawData, dataDir)
        fileList = [f for f in os.listdir(fullDataDir) if os.path.isfile(os.path.join(fullDataDir, f))]

        # Find the obsId and aFilter
        obsIdList = []
        for fileName in fileList:
            m = re.match(r"\S*_(\d*)_f(\d)_\S*", fileName)
            if (m is not None):
                obsIdList.append(m.groups())

        # Get the unique list
        obsIdList = list(set(obsIdList))
        if (len(obsIdList) > 1):
            raise RuntimeError("There are more than one unique ObdId and filter in directory.")
        data = obsIdList[0]

        # Import to butler
        phosimFilterID = {"0": "u", "1": "g", "2": "r", "3": "i", "4": "z", "5": "y"}
        obsId = int(data[0])
        aFilter = phosimFilterID[data[1]]
        self.dataCollector.importPhoSimDataToButler(dataDir, obsId=obsId, aFilter=aFilter, atype=atype, 
                                                    overwrite=overwrite)

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

    def doISR(self, visit, sensorName, snap=0, fakeDatasetType="eimage", 
                outputDatasetType="postISRCCD"):
        """
        
        Do the instrument signature removal (ISR).
        
        Arguments:
            visit {[int]} -- Visit time.
            sensorName {[str]} -- Sensor name. (e.g. "R:2,2 S:1,1")
        
        Keyword Arguments:
            snap {int} -- Snap time (0 or 1) means first/ second exposure. (default: {0})
            fakeDatasetType {[str]} -- Use this type of image supported by lsst camera mapper 
                                        to simulate the post-ISR image. (default: {"eimage"})
            outputDatasetType {[str]} -- Output data type supported by lsst camera mapper. 
                                        (default: {"postISRCCD"})

        Returns:
            [ExposureU] -- Exposure image after ISR.
        """

        # Use the regular expression to analyze the input name
        m = re.match(r"R:(\d,\d) S:(\d,\d)(?:,([A,B]))?$", sensorName)
        if (m is not None):
            raft, sensor = m.groups()[0:2]

            # Do the ISR
            self.isrWrapper.doISR(visit, snap, raft, sensor, channel=None, 
                        fakeDatasetType=fakeDatasetType, outputDatasetType=outputDatasetType)
        else:
            raise RuntimeError("Sensor name: '%s' is not allowed." % sensorName)

if __name__ == "__main__":
    
    # Instintiate the components
    sourSelc = SourceSelector()
    dataCollector = WFDataCollector()
    isrWrapper = EimgIsrWrapper()

    # Configurate the source selector
    cameraType = "comcam"
    dbType = "LocalDb"
    aFilter = "g"
    cameraMJD = 59580.0

    sourSelc.configSelector(cameraType=cameraType, dbType=dbType, aFilter=aFilter, 
                            cameraMJD=cameraMJD)

    # Set the criteria of neighboring stars
    starRadiusInPixel = 63
    spacingCoefficient = 2.5
    sourSelc.configNbrCriteria(starRadiusInPixel, spacingCoefficient)

    # Configurate the wfs data collector
    pathOfRawData = "../test/phosimOutput"
    destinationPath = "../test"
    butlerInputs = "../test"
    butlerOutputs = "../test"
    regisAdress = "../test/registry.sqlite3"
    dataCollector.config(pathOfRawData=pathOfRawData, destinationPath=destinationPath, 
                dbAdress=regisAdress, butlerInputs=butlerInputs, butlerOutputs=butlerOutputs)

    # Configurate the ISR wrapper
    isrWrapper.configWrapper(inputs=butlerInputs, outputs=butlerOutputs)

    # Initiate the WEP Controller
    wepCntlr = WEPController()
    wepCntlr.config(sourSelc=sourSelc, dataCollector=dataCollector, isrWrapper=isrWrapper)

    # Set the database address
    dbAdress = "../test/bsc.db3"

    # Do the query
    pointing = (0,0)
    cameraRotation = 0.0
    skyInfoFilePath = "../test/phosimOutput/realComCam/output/skyComCamInfo.txt"

    neighborStarMap, starMap, wavefrontSensors = wepCntlr.getTargetStarByFile(dbAdress, skyInfoFilePath, 
                                        pointing, cameraRotation, orientation="all", tableName="TempTable")

    # Import the PhoSim simulated image
    dataDirList = ["realComCam/output/Extra", "realComCam/output/Intra"]
    # for dataDir in dataDirList:
    #     wepCntlr.importPhoSimDataToButler(dataDir, atype="raw", overwrite=False)

    # Do the ISR
    extraObsId = 9007000
    intraObsId = 9007001
    obsIdList = [extraObsId, intraObsId]
    sensorNameList = list(starMap.keys())
    for obsId in obsIdList:
        for sensorName in sensorNameList:
            wepCntlr.doISR(obsId, sensorName)








    # # Import the PhoSim simulated image
    # extraObsId = 9007000
    # intraObsId = 9007001
    # obsIdList = [extraObsId, intraObsId]
    # aFilter = "g"

    # # Set the sensor information
    # snap = 0
    # raft = "2,2"
    # sensor = "1,1"

    # # Do the fake ISR
    # wepCntlr.configIsrWrapper(inputs=butlerInputs, outputs=butlerOutputs)

    # for ii in range(2):
    #     wepCntlr.doISR(obsIdList[ii], snap, raft, sensor, fakeDatasetType="eimage", 
    #                     outputDatasetType="postISRCCD")

    # # Get the image data
    # intraImg, extraImg = wepCntlr.getDefocalImg(snap, raft, sensor, intraObsId, extraObsId, 
    #                                         datasetType="postISRCCD")