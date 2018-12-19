import os
import unittest

from lsst.ts.wep.LocalDatabaseDecorator import LocalDatabaseDecorator
from lsst.ts.wep.Utility import getModulePath


class TestLocalDatabaseDecorator(unittest.TestCase):
    """Test the local database decorator class."""

    def setUp(self):

        # Get the path of module
        self.modulePath = getModulePath()

        # Address of local database
        dbAdress = os.path.join(self.modulePath, "tests", "testData",
                                "bsc.db3")

        self.db = LocalDatabaseDecorator()
        self.db.connect(dbAdress)

    def tearDown(self):
        self.db.disconnect()

    def testFunctions(self):
        
        aFilter = "g"
        tableName = "TempTable"
        self.assertFalse(self.db.checkTableInDb(tableName))

        self.db.createTable(aFilter, tableName)
        self.assertTrue(self.db.checkTableInDb(tableName))

        try:
            self.db.createTable(aFilter, tableName)
        except Exception as RuntimeError:
            pass

        # Sky data file path
        skyFilePath = os.path.join(self.modulePath, "tests", "testData",
                                   "skyComCamInfo.txt")
        self.db.insertDataByFile(aFilter, tableName, skyFilePath)

        self.db.deleteTable(tableName)
        self.assertFalse(self.db.checkTableInDb(tableName))


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
