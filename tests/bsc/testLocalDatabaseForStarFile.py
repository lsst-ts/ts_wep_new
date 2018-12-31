import os
import unittest

from lsst.ts.wep.bsc.LocalDatabaseForStarFile import LocalDatabaseForStarFile
from lsst.ts.wep.Utility import getModulePath, FilterType


class TestLocalDatabaseForStarFile(unittest.TestCase):
    """Test the local database for star file class."""

    def setUp(self):

        self.filterType = FilterType.G
        self.db = LocalDatabaseForStarFile()

        self.modulePath = getModulePath()
        dbAdress = os.path.join(self.modulePath, "tests", "testData",
                                "bsc.db3")
        self.db.connect(dbAdress)

    def tearDown(self):

        self.db.deleteTable(self.filterType)
        self.db.disconnect()

    def testTableIsInDb(self):

        self.assertFalse(self.db._tableIsInDb("StarTableG"))
        self.assertTrue(self.db._tableIsInDb("BrightStarCatalogU"))

    def testCreateTable(self):

        self._createTable()
        self.assertTrue(self.db._tableIsInDb("StarTableG"))

    def _createTable(self):
        self.db.createTable(self.filterType)

    def testCreateTableIfTableExist(self):

        self._createTable()
        self.assertRaises(ValueError, self.db.createTable, FilterType.G)

    def testInsertDataByFile(self):

        self._createTable()
        idAll = self.db.getAllId(self.filterType)
        self.assertEqual(len(idAll), 0)

        skyFilePath = os.path.join(self.modulePath, "tests", "testData",
                                   "skyComCamInfo.txt")
        self.db.insertDataByFile(skyFilePath, self.filterType)
        idAll = self.db.getAllId(self.filterType)

        self.assertEqual(len(idAll), 4)

    def testDeleteTable(self):

        self._createTable()
        self.assertTrue(self.db._tableIsInDb("StarTableG"))

        self.db.deleteTable(self.filterType)
        self.assertFalse(self.db._tableIsInDb("StarTableG"))


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
