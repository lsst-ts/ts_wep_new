import os
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
        
        self.dataCollector = WFDataCollector()
        self.isrWrapper = EimgIsrWrapper()
        self.sourSelc = SourceSelector()
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

    def setFilter(self, atype):
        """
        
        Set the active filter type.
        
        Arguments:
            atype {[str]} -- Filter type ("u", "g", "r", "i", "z", "y").
        """

        self.sourSelc.filter.setFilter(atype)

    def connectDb(self, *kwargs):
        """
        
        Connect the database.
        
        Arguments:
            *kwargs {[string]} -- Information to connect to the database.
        """

        self.sourSelc.connect(*kwargs)

    def disconnectDb(self):
        """
        
        Disconnect the database.
        """

        self.sourSelc.disconnect()

    def getTargetStarByFile(self, dbAdress, skyInfoFilePath, pointing, cameraRotation, 
                            orientation=None, offset=0, tableName="TempTable"):
        """
        
        Get the target stars by querying the file.
        
        Arguments:
            dbAdress {[str]} -- Local database address.
            skyInfoFilePath {[str]} -- File path of sky information.
            pointing {[tuple]} -- Camera boresight (RA, Decl) in degree.
            cameraRotation {[float]} -- Camera rotation angle in degree.
        
        Keyword Arguments:
            orientation {[str]} -- Orientation of wavefront sensor(s) on camera. (default: {None})
            offset {[float]} -- Add offset in pixel to sensor dimension for judging stars on detector or not. 
                                offset=0 for normal use. 
                                offset=maxDistance to generate the local database. (default: {0})
            tableName {[str]} -- Table name. (default: {None})
        
        Returns:
            neighborStarMap {[list]} -- Information of neighboring stars and candidate stars with 
                                       the name of sensor as a list.
            starMap {[list]} -- Information of stars with the name of sensor as a list.
            wavefrontSensors {[list]} -- Corners of sensor with the name of sensor as a list.
        """

        # Create the table
        localDb = LocalDatabaseDecorator()
        localDb.connect(dbAdress)

        # Insert the sky data
        aFilter = self.sourSelc.filter.getFilter()
        localDb.createTable(aFilter, tableName)
        localDb.insertDataByFile(aFilter, tableName, skyInfoFilePath)
        localDb.disconnect()
        
        # Do the query and analysis
        self.connectDb(dbAdress)
        neighborStarMap, starMap, wavefrontSensors = self.getTargetStar(pointing, cameraRotation, 
                                        orientation=orientation, offset=offset, tableName=tableName)
        self.disconnectDb()

        # Delete the table
        localDb.connect(dbAdress)
        localDb.deleteTable(tableName)
        localDb.disconnect()

        return neighborStarMap, starMap, wavefrontSensors

    def getTargetStar(self, pointing, cameraRotation, orientation=None, offset=0, tableName=None):
        """
        
        Get the target stars by querying the database.
        
        Arguments:
            pointing {[tuple]} -- Camera boresight (RA, Decl) in degree.
            cameraRotation {[float]} -- Camera rotation angle in degree.
        
        Keyword Arguments:
            orientation {[str]} -- Orientation of wavefront sensor(s) on camera. (default: {None})
            offset {[float]} -- Add offset in pixel to sensor dimension for judging stars on detector or not. 
                                offset=0 for normal use. 
                                offset=maxDistance to generate the local database. (default: {0})
            tableName {[str]} -- Table name. (default: {None})
        
        Returns:
            neighborStarMap {[list]} -- Information of neighboring stars and candidate stars with 
                                       the name of sensor as a list.
            starMap {[list]} -- Information of stars with the name of sensor as a list.
            wavefrontSensors {[list]} -- Corners of sensor with the name of sensor as a list.
        """

        neighborStarMap, starMap, wavefrontSensors = self.sourSelc.getTargetStar(pointing, 
                    cameraRotation, orientation=orientation, offset=offset, tableName=tableName)

        return neighborStarMap, starMap, wavefrontSensors

    def configNeiborStarCriteria(self, starRadiusInPixel, spacingCoefficient, maxNeighboringStar=99):
        """
        
        Set the configuration to decide the scientific target.
        
        Arguments:
            starRadiusInPixel {[float]} -- Diameter of star. For the defocus = 1.5 mm, the star's 
                                            radius is 63 pixel.
            spacingCoefficient {[float]} -- Maximum distance in units of radius one donut must be 
                                            considered as a neighbor.
        
        Keyword Arguments:
            maxNeighboringStar {[int]} -- Maximum number of neighboring stars. (default: {99})
        """

        self.sourSelc.config(starRadiusInPixel, spacingCoefficient, 
                                maxNeighboringStar=maxNeighboringStar)

    def configSourSelc(self, cameraType, dbType=None, cameraMJD=59580.0):
        """
        
        Do the configuration of source selector.
        
        Arguments:
            cameraType {[str]} -- Type of camera ("lsst" or "comcam").
        
        Keyword Arguments:
            dbType {[str]} -- Type of database ("UWdb" or "LocalDb"). (default: {None}) 
            cameraMJD {float} -- Camera MJD. (default: {59580.0})
        """
        
        self.sourSelc.configSelector(cameraType=cameraType, dbType=dbType, cameraMJD=cameraMJD)

    def configIsrWrapper(self, inputs=None, outputs=None):
        """
        
        Do the configuration of instrument signature remover (ISR) wrapper.
        
        Keyword Arguments:
            inputs {[RepositoryArg or string]} -- Can be a single item or a list. Provides arguments 
                    to load an existing repository (or repositories). String is assumed to be a URI 
                    and is used as the cfgRoot (URI to the location of the cfg file). (Local file system 
                    URI does not have to start with 'file://' and in this way can be a relative path).
                    (default: {None})
            outputs {[RepositoryArg or string]} -- Can be a single item or a list. Provides arguments to 
                    load one or more existing repositories or create new ones. String is assumed to be a 
                    URI and as used as the repository root.
        """

        self.isrWrapper.configWrapper(inputs=inputs, outputs=outputs)

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

    def doISR(self, visit, snap, raft, sensor, channel=None, fakeDatasetType="eimage", 
                outputDatasetType="postISRCCD"):
        """
        
        Do the instrument signature removal (ISR).
        
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
            [ExposureU] -- Exposure image after ISR.
        """

        return self.isrWrapper.doISR(visit, snap, raft, sensor, channel=channel, 
                    fakeDatasetType=fakeDatasetType, outputDatasetType=outputDatasetType)

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
        wepCntlr.importPhoSimDataToButler(dataDirList[ii], obsId=obsIdList[ii], aFilter=aFilter, 
                                            atype=atype, overwrite=False)

    # Set the sensor information
    snap = 0
    raft = "2,2"
    sensor = "1,1"

    # Do the fake ISR
    wepCntlr.configIsrWrapper(inputs=butlerInputs, outputs=butlerOutputs)

    for ii in range(2):
        wepCntlr.doISR(obsIdList[ii], snap, raft, sensor, fakeDatasetType="eimage", 
                        outputDatasetType="postISRCCD")

    # Get the image data
    intraImg, extraImg = wepCntlr.getDefocalImg(snap, raft, sensor, intraObsId, extraObsId, 
                                            datasetType="postISRCCD")

    # Do the source selector
    cameraMJD = 59580.0
    cameraType = "comcam"
    wepCntlr.configSourSelc(cameraType, dbType="LocalDb", cameraMJD=cameraMJD)

    # Set the database address
    dbAdress = "../test/bsc.db3"

    # Do the query
    pointing = (0,0)
    cameraRotation = 0.0
    skyInfoFilePath = "../test/phosimOutput/realComCam/output/skyComCamInfo.txt"

    starRadiusInPixel = 63
    spacingCoefficient = 2.5

    wepCntlr.setFilter(aFilter)
    wepCntlr.configNeiborStarCriteria(starRadiusInPixel, spacingCoefficient)

    # wepCntlr.connectDb(dbAdress)
    # neighborStarMap, starMap, wavefrontSensors = wepCntlr.getTargetStar(pointing, cameraRotation, orientation="all")
    # wepCntlr.disconnectDb()

    neighborStarMap, starMap, wavefrontSensors = wepCntlr.getTargetStarByFile(dbAdress, skyInfoFilePath, pointing, 
                                                                                cameraRotation, orientation="all")