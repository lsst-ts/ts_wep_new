import os
import unittest

from lsst.ts.wep.bsc.StarData import StarData
from lsst.ts.wep.bsc.LocalDatabase import LocalDatabase
from lsst.ts.wep.Utility import getModulePath, FilterType


class TestLocalDatabase(unittest.TestCase):
    """Test the LocalDatabase class."""

    def setUp(self):

        # Set the database address
        modulePath = getModulePath()
        dbAdress = os.path.join(modulePath, "tests", "testData", "bsc.db3")

        # Set up local database
        self.localDatabase = LocalDatabase()
        self.localDatabase.connect(dbAdress)

        # Set the filter
        self.filterType = FilterType.G

        # Set up neighboring star map
        stars = StarData([123, 456, 789], [0.1, 0.2, 0.3], [2.1, 2.2, 2.3],
                         [2.0, 3.0, 4.0], [2.0, 3.0, 4.0], [2.0, 3.0, 4.0],
                         [2.0, 3.0, 4.0], [2.0, 3.0, 4.0], [2.0, 3.0, 4.0])
        stars.setRaInPixel(stars.getRA() * 10)
        stars.setDeclInPixel(stars.getDecl() * 10)
        self.neighboringStar = stars.getNeighboringStar(
            [0], 3, self.filterType, 99)

    def tearDown(self):

        allId = self.localDatabase.getAllId(self.filterType)
        self.localDatabase.deleteData(self.filterType, allId)

        self.localDatabase.disconnect()

    def testQuery(self):

        self._insertData()
        stars = self.localDatabase.query(self.filterType, [0, 2],
                                         [0, 2.4], [0.4, 2],
                                         [0.4, 2.4])

        self.assertEqual(stars.getId().tolist(), [123, 456, 789])
        self.assertEqual(stars.getRA().tolist(), [0.1, 0.2, 0.3])
        self.assertEqual(stars.getDecl().tolist(), [2.1, 2.2, 2.3])

    def testQueryWithoutStarAndCrossRa0(self):

        stars = self._queryCrossRa0()
        self.assertEqual(len(stars.getId()), 0)

    def testQueryWithStarAndCrossRa0(self):

        self._insertData()
        listID = self._getListId()
        self.localDatabase.updateData(self.filterType, listID,
                                      ["ra", "ra", "ra"],
                                      [0.005, 359.98, 359.999])
        self.localDatabase.updateData(self.filterType, listID,
                                      ["decl", "decl", "decl"],
                                      [-1.5, -1.5, -1.5])

        stars = self._queryCrossRa0()

        self.assertEqual(len(stars.getId()), 2)
        self.assertEqual(stars.getId().tolist(), [123, 789])
        self.assertEqual(stars.getRA().tolist(), [0.005, 359.999])
        self.assertEqual(stars.getDecl().tolist(), [-1.5, -1.5])

    def _queryCrossRa0(self):

        stars = self.localDatabase.query(self.filterType, [0.01, -1],
                                         [359.99, -2], [0.01, -1],
                                         [359.99, -2])
        return stars

    def testSearchSimobjdID(self):

        starData = self.localDatabase.searchSimobjdID(FilterType.U,
                                                      [54408946])
        self.assertEqual(len(starData), 1)
        self.assertEqual(starData[0][0], 2)
        self.assertEqual(starData[0][1], 359.732296)
        self.assertEqual(starData[0][2], 63.053469)

    def testSearchRaDecl(self):

        starId = self.localDatabase.searchRaDecl(FilterType.U, 359.732296,
                                                 63.053469)
        self.assertEqual(starId[0], 2)

    def testInsertData(self):

        self._insertData()
        starId = self.localDatabase.searchRaDecl(self.filterType, 0.1, 2.1)
        self.assertEqual(len(starId), 1)

    def testInsertDataRepeatly(self):

        for ii in range(5):
            self._insertData()

        allId = self.localDatabase.getAllId(self.filterType)
        self.assertEqual(len(allId), 3)

    def _insertData(self):

        self.localDatabase.insertData(self.filterType, self.neighboringStar)

    def testUpdateData(self):

        self._insertData()
        listID = self._getListId()
        self.localDatabase.updateData(self.filterType, listID,
                                      ["ra", "ra", "ra"], [1.0, 2.0, 3.0])

        oldDataId = self.localDatabase.searchRaDecl(self.filterType,
                                                    0.2, 2.2)
        newDataId = self.localDatabase.searchRaDecl(self.filterType,
                                                    2.0, 2.2)
        self.assertEqual(len(oldDataId), 0)
        self.assertEqual(len(newDataId), 1)

    def _getListId(self):

        starData = self.localDatabase.searchSimobjdID(self.filterType,
                                                      [123, 456, 789])
        listID = [starData[0][0], starData[1][0], starData[2][0]]

        return listID

    def testDeleteData(self):

        self._insertData()
        listID = self._getListId()
        self.localDatabase.deleteData(self.filterType, listID)

        allId = self.localDatabase.getAllId(self.filterType)
        self.assertEqual(len(allId), 0)

    def testGetAllIdWithData(self):

        allId = self.localDatabase.getAllId(FilterType.U)
        self.assertEqual(len(allId), 7155)

    def testGetAllIdWithoutData(self):

        allId = self.localDatabase.getAllId(self.filterType)
        self.assertEqual(len(allId), 0)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
