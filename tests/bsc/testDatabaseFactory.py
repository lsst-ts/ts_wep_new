import unittest

from lsst.ts.wep.bsc.DatabaseFactory import DatabaseFactory
from lsst.ts.wep.bsc.LocalDatabase import LocalDatabase
from lsst.ts.wep.bsc.LocalDatabaseForStarFile import LocalDatabaseForStarFile
from lsst.ts.wep.Utility import BscDbType


class TestDatabaseFactory(unittest.TestCase):
    """Test the DatabaseFactory class."""

    def testCreateDb(self):

        db = DatabaseFactory.createDb(BscDbType.LocalDb)
        self.assertTrue(isinstance(db, LocalDatabase))

        dbForFile = DatabaseFactory.createDb(BscDbType.LocalDbForStarFile)
        self.assertTrue(isinstance(dbForFile, LocalDatabaseForStarFile))

        self.assertRaises(ValueError, DatabaseFactory.createDb, "wrongType")


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
