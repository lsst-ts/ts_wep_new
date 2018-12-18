import os
import unittest

from lsst.ts.wep.isr.LocalDatabase import LocalDatabase
from lsst.ts.wep.Utility import getModulePath


class TestLocalDatabase(unittest.TestCase):
    """Test the LocalDatabase class."""

    def setUp(self):

        # Get the path of module
        modulePath = getModulePath()

        # Local database setting
        dbAdress = os.path.join(modulePath, "tests", "testData",
                                "registry.sqlite3")
        self.dbAdress = dbAdress

    def testLocalDatabase(self):

        # Connect the local database
        localDatabase = LocalDatabase()
        localDatabase.connect(self.dbAdress)

        # Parameters from Phosim
        obsId = 99999998
        snap = 0
        raft = "2,2"
        sensor = "1,1"
        aFilter = "r"

        # Pretend to have a visit time (tai: International Atomic Time)
        taiObs = "1994-07-19T06:49:51.543999744"

        # Pretend to have the skyTile
        skyTile = 30001

        # Exposure time is 15 sec in single snap
        expTime = 15.0

        # Test to query a certain visit
        item = localDatabase.query("raw", visit=obsId, snap=snap, raft=raft,
                                   sensor=sensor, aFilter=aFilter)
        self.assertEqual(item, [])

        # Test to insert the data
        localDatabase.insertData(obsId, snap, raft, sensor, aFilter, taiObs,
                                 skyTile, expTime)

        # Query the inserted data
        item = localDatabase.query("raw", visit=obsId, snap=snap, raft=raft,
                                   sensor=sensor, aFilter=aFilter)
        self.assertNotEqual(item, [])

        # Test to delete the data
        localDatabase.deleteData(obsId, snap, raft, sensor, aFilter)

        # Query the inserted data
        item = localDatabase.query("raw", visit=obsId, snap=snap, raft=raft,
                                   sensor=sensor, aFilter=aFilter)
        self.assertEqual(item, [])

        # Disconnect from the local database
        localDatabase.disconnect()


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
