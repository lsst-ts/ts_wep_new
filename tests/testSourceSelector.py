import os
import numpy as np
import unittest

from lsst.ts.wep.SourceSelector import SourceSelector, calcPixPos
from lsst.ts.wep.Utility import getModulePath


class TestSourceSelector(unittest.TestCase):
    """Test the source selector class."""

    def setUp(self):

        # Get the path of module
        self.modulePath = getModulePath()

        # Camera type: "lsst" or "comcam"
        cameraType = "comcam"

        # Active filter type
        aFilterType = "r"

        # Address of local database
        dbAdress = os.path.join(self.modulePath, "tests", "testData",
                                "bsc.db3")

        # Remote database setting
        databaseHost = "localhost:51433"
        databaseUser = "LSST-2"
        databasePassword = "L$$TUser"
        databaseName = "LSSTCATSIM"

        # Set the database
        self.remoteDb = SourceSelector()
        self.localDb = SourceSelector()

        self.remoteDb.configSelector(cameraType=cameraType, dbType="UWdb",
                                     aFilter=aFilterType)
        self.localDb.configSelector(cameraType=cameraType, dbType="LocalDb",
                                    aFilter=aFilterType)

        # Remote database infomation
        remoteDbInfo = [databaseHost, databaseUser, databasePassword,
                        databaseName]

        # Connect to database
        self.remoteDb.connect(*remoteDbInfo)
        self.localDb.connect(dbAdress)

    def tearDown(self):

        # Disconnect database
        self.remoteDb.disconnect()
        self.localDb.disconnect()

    def testFunctions(self):

        # Boresight (RA, Dec) (unit: degree) (0 <= RA <= 360, -90 <= Dec <= 90)
        pointing = (20.0, 30.0)

        # Camera rotation
        cameraRotation = 0.0

        # Camera orientation for ComCam ("center" or "corner" or "all")
        # Camera orientation for LSSTcam ("corner" or "all")
        orientation = "center"

        # Maximum distance in units of radius one donut must be considered as
        # a neighbor.
        spacingCoefficient = 2.5

        # For the defocus = 1.5 mm, the star's radius is 63 pixel.
        starRadiusInPixel = 63

        # Set the configuration to select the scientific target
        self.remoteDb.configNbrCriteria(starRadiusInPixel, spacingCoefficient,
                                        maxNeighboringStar=99)
        self.localDb.configNbrCriteria(starRadiusInPixel, spacingCoefficient,
                                       maxNeighboringStar=99)

        # Set the active filter
        # self.remoteDb.setFilter(self.aFilterType)
        # self.localDb.setFilter(self.aFilterType)

        # Test to get the active filter
        self.assertEqual(self.localDb.getFilter(), "r")

        # Test to get the standard deviation 
        self.assertEqual(self.localDb.getStddevSplit(), 20.0)

        # Get the scientific target by querying the remote database
        neighborStarMap, starMap, wavefrontSensors = \
            self.remoteDb.getTargetStar(pointing, cameraRotation,
                                        orientation=orientation)

        # Test to get at least one star
        allStars = starMap["R:2,2 S:1,1"]
        self.assertTrue(len(allStars.SimobjID)>=1)

        # Get the scientific target by querying the local database
        neighborStarMapLocal, starMapLocal, wavefrontSensorsLocal = \
            self.localDb.getTargetStar(pointing, cameraRotation,
                                       orientation=orientation)

        # Test the get the empty star map
        allStarsLocal = starMapLocal["R:2,2 S:1,1"]
        self.assertEqual(allStarsLocal.SimobjID, [])

        # Insert the neighboring star map into the database
        self.localDb.insertToBSC(neighborStarMap)

        # Query the local database again
        neighborStarMapLocal, starMapLocal, wavefrontSensorsLocal = \
            self.localDb.getTargetStar(pointing, cameraRotation,
                                       orientation=orientation)

        # Test to get all neighboring stars
        allNeighborStarLocal = neighborStarMapLocal["R:2,2 S:1,1"]
        allNeighborStar = neighborStarMap["R:2,2 S:1,1"]
        self.assertEqual(len(allNeighborStarLocal.SimobjID),
                         len(allNeighborStar.SimobjID))

        # Test to trim the margin
        self.remoteDb.trimMargin(neighborStarMap, 1000)

        # Test to search the id of star based on (ra, decl)
        searchStarId = self.localDb.searchRaDecl(20.088157, 29.983533)

        # Test to update the value
        self.localDb.updateBSC([searchStarId[0][0], searchStarId[0][0]],
                               ["ra", "decl"], [200, 200])
        newSearchStarId = self.localDb.searchRaDecl(200, 200)
        self.assertEqual(searchStarId, newSearchStarId)

        # Delete all data in local database
        allStarList = np.arange(1,len(allNeighborStarLocal.RaDecl)+1)
        self.localDb.db.deleteData(self.localDb.getFilter(), allStarList.tolist())

    def testCoorFun(self):
        fitsFilePath = os.path.join(self.modulePath, "tests", "testData",
                                    "eimage", "v99999999-fr", "E000", "R22",
                                    "eimage_99999999_R22_S11_E000.fits.gz")
        raList = [0]
        decList = [0]
        xPosList, yPosList = calcPixPos(fitsFilePath, raList, decList)
        self.assertEqual((xPosList[0], yPosList[0]), (2000, 2036))


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
