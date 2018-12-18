from decimal import Decimal
import unittest

from lsst.ts.wep.bsc.StarData import StarData
from lsst.ts.wep.bsc.BrightStarDatabase import BrightStarDatabase


class TestBrightStarDatabase(unittest.TestCase):
    """Test the BrightStarDatabase class."""

    # UW database setting
    databaseHost = "localhost:51433"
    databaseUser = "LSST-2"
    databasePassword = "L$$TUser"
    databaseName = "LSSTCATSIM"
 
    # Camera Filter
    cameraFilter = "g"

    # Neighboring stars
    neighboringStar = None

    def setUp(self):
        # Set up UW database
        self.database = BrightStarDatabase()
        self.database.connect(self.databaseHost, self.databaseUser, 
                              self.databasePassword, self.databaseName)
        
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
        self.database.disconnect()
        
    def testFoobar(self):
        stars = self.database.query("bright_stars", self.database.FilterU,
                                    [0.01, -1], [359.99, -2], [0.01, -1],
                                    [359.99, -2])
        
        self.assertEqual(len(stars.LSSTMagU), 36)
        self.assertEqual(max(stars.LSSTMagU), 30.8771)
        self.assertEqual(min(stars.LSSTMagU), 14.74348)

    def testUQuery(self):
        stars = self.database.query("bright_stars", self.database.FilterU,
                                    [75.998622, -1], [75.998622, -2],
                                    [75.998985, -1], [75.998985, -2])

        self.assertEqual(len(stars.SimobjID), 3)
        self.assertEqual(len(stars.RA), 3)
        self.assertEqual(len(stars.Decl), 3)
        self.assertEqual(len(stars.LSSTMagU), 3)

        self.assertEqual(stars.SimobjID[0], Decimal("375921344"))
        self.assertEqual(stars.SimobjID[1], Decimal("377778621"))
        self.assertEqual(stars.SimobjID[2], Decimal("378267883"))

        self.assertEqual(stars.RA[0], 75.998623)
        self.assertEqual(stars.RA[1], 75.998787)
        self.assertEqual(stars.RA[2], 75.998985)

        self.assertEqual(stars.Decl[0], -1.526383)
        self.assertEqual(stars.Decl[1], -1.15006)
        self.assertEqual(stars.Decl[2], -1.007737)

        self.assertEqual(stars.LSSTMagU[0], 17.4837)
        self.assertEqual(stars.LSSTMagU[1], 24.45101)
        self.assertEqual(stars.LSSTMagU[2], 15.53796)

        simobjid = self.database.searchRaDecl("bright_stars", 359.432736,
                                              -80.315527)
        self.assertEqual(simobjid[0][0], Decimal("5760469"))


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
