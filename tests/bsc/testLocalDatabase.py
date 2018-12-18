import os
import unittest

from lsst.ts.wep.bsc.StarData import StarData
from lsst.ts.wep.bsc.LocalDatabase import LocalDatabase
from lsst.ts.wep.Utility import getModulePath


class TestLocalDatabase(unittest.TestCase):
    """Test the LocalDatabase class."""

    # Camera Filter
    cameraFilter = "g"

    # Neighboring stars
    neighboringStar = None

    def setUp(self):

        # Get the path of module
        modulePath = getModulePath()

        # Local database setting
        dbAdress = os.path.join(modulePath, "tests", "testData", "bsc.db3")
        
        # Set up local database
        self.localDatabase = LocalDatabase()
        self.localDatabase.connect(dbAdress)

        # Set up neighboring star map
        stars = StarData([123, 456, 789], [0.1, 0.2, 0.3], [2.1, 2.2, 2.3],
                         [2.0, 3.0, 4.0], [2.0, 3.0, 4.0], [2.0, 3.0, 4.0],
                         [2.0, 3.0, 4.0], [2.0, 3.0, 4.0], [2.0, 3.0, 4.0])
        stars.populateRAData([value*10 for value in stars.RA])
        stars.populateDeclData([value*10 for value in stars.Decl])
        self.neighboringStar = stars.getNeighboringStar(
                                            [0], 3, self.cameraFilter, 99)

    def tearDown(self):
        # Disconnect database
        self.localDatabase.disconnect()
        
    def testLocalDatabase(self):

        localTableName = "BrightStarCatalog" + self.cameraFilter.upper()

        # Test to add data
        self.localDatabase.insertData(self.cameraFilter, self.neighboringStar)
        id = self.localDatabase.searchRaDecl(self.cameraFilter, 0.1, 2.1)
        self.assertEqual(len(id), 1)

        # Test to add the repeated data
        numId = self.localDatabase.getAllId(self.cameraFilter)
        self.localDatabase.insertData(self.cameraFilter, self.neighboringStar)
        newNumId = self.localDatabase.getAllId(self.cameraFilter)
        self.assertEqual(len(numId), len(newNumId))

        # Test to search simobjdID
        simobjid = self.localDatabase.searchSimobjdID(self.cameraFilter,
                                                      [123, 456, 789])
        self.assertEqual(len(simobjid), 3)

        # Test to update data
        self.localDatabase.updateData(
                            self.cameraFilter,
                            [simobjid[0][0], simobjid[1][0], simobjid[2][0]], 
                            ["ra", "ra", "ra"], [1.0, 2.0, 3.0])
        oldDataId = self.localDatabase.searchRaDecl(self.cameraFilter,
                                                    0.2, 2.2)
        newDataId = self.localDatabase.searchRaDecl(self.cameraFilter,
                                                    2.0, 2.2)
        self.assertEqual(len(oldDataId), 0)
        self.assertEqual(len(newDataId), 1)

        # Test to delete data
        self.localDatabase.deleteData(
            self.cameraFilter, [simobjid[0][0], simobjid[1][0], simobjid[2][0]])
        newNumId = self.localDatabase.getAllId(self.cameraFilter)
        self.assertNotEqual(numId, newNumId)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
