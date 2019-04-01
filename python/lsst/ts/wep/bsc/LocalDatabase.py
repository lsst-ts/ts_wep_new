import sqlite3
import numpy as np

from lsst.ts.wep.bsc.DefaultDatabase import DefaultDatabase
from lsst.ts.wep.bsc.StarData import StarData
from lsst.ts.wep.Utility import FilterType


class LocalDatabase(DefaultDatabase):

    PRE_TABLE_NAME = "BrightStarCatalog"

    def connect(self, dbAdress):
        """Connects database based on the local path.

        Parameters
        ----------
        dbAdress : str
            Path of local sqlite3 database.
        """

        self.connection = sqlite3.connect(dbAdress)
        self.cursor = self.connection.cursor()

    def _queryTable(self, filterType, top, bottom, left, right):
        """Queries the database for stars within an area.

        Parameters
        ----------
        filterType : FilterType
            Filter type.
        top : float
            The top edge of the box (Decl).
        bottom : float
            The bottom edge of the box (Decl).
        left : float
            The left edge of the box (RA).
        right : float
            The right edge of the box (RA).

        Returns
        ----------
        StarData
            Star information.
        """

        # Do the query
        tableName = self._getTableName(filterType)
        command = "SELECT simobjid, ra, decl, " + \
                  filterType.name.lower() + "mag" + \
                  " FROM " + tableName + \
                  " WHERE decl <= %f AND decl >= %f AND ra >= %f AND ra <= %f"
        query = command % (top, bottom, left, right)
        self.cursor.execute(query)

        # Collect the data
        simobjid = []
        ra = []
        decl = []
        lsstMagU = []
        lsstMagG = []
        lsstMagR = []
        lsstMagI = []
        lsstMagZ = []
        lsstMagY = []
        for item in self.cursor.fetchall():

            # It is noted that the data type of simobjid is big interger
            # in UW database
            simobjid.append(item[0])
            ra.append(item[1])
            decl.append(item[2])

            if (filterType == FilterType.U):
                lsstMagU.append(item[3])

            elif (filterType == FilterType.G):
                lsstMagG.append(item[3])

            elif (filterType == FilterType.R):
                lsstMagR.append(item[3])

            elif (filterType == FilterType.I):
                lsstMagI.append(item[3])

            elif (filterType == FilterType.Z):
                lsstMagZ.append(item[3])

            elif (filterType == FilterType.Y):
                lsstMagY.append(item[3])

        return StarData(simobjid, ra, decl, lsstMagU, lsstMagG, lsstMagR,
                        lsstMagI, lsstMagZ, lsstMagY)

    def _getTableName(self, filterType):
        """Get the table name.

        Parameters
        ----------
        filterType : FilterType
            Filter type.

        Returns
        -------
        str
            Table name.
        """

        return self.PRE_TABLE_NAME + filterType.name

    def searchSimobjdID(self, filterType, listID):
        """Search the data based on the simobjid.

        Parameters
        ----------
        filterType : FilterType
            Filter type.
        listID : list[int]
            Simobjid list to search.

        Returns
        -------
        list[tuple]
            Results [(id_1, ra_1, decl_1, ), (id_2, ra_2, decl_2, ), ...] of
            search
        """

        # Search the simobjid data
        tableName = self._getTableName(filterType)
        command = "SELECT id, ra, decl From " + tableName + \
                  " WHERE simobjid in" + \
                  " (" + ', '.join(str(ID) for ID in listID) + ")"
        self.cursor.execute(command)

        return self.cursor.fetchall()

    def searchRaDecl(self, filterType, ra, decl):
        """Search the star Id based on ra, decl.

        Parameters
        ----------
        filterType : FilterType
            Filter type.
        ra : float
            Star right ascension in degree.
        decl : float
            Star declination in degree.

        Returns
        -------
        list
            Star Id list in local database.
        """

        # Compare ra and decl to see the existance of star in database
        tableName = self._getTableName(filterType)
        command = "SELECT id FROM " + tableName + \
                  " WHERE ra = %f AND decl = %f"
        query = command % (ra, decl)
        self.cursor.execute(query)

        listIdByQuery = self.cursor.fetchall()

        return self._changeQueriedIdToList(listIdByQuery)

    def _changeQueriedIdToList(self, listIdByQuery):
        """Change the queried Id to the list data type.

        Parameters
        ----------
        listIdByQuery : list[tuple]
            Queried Id list.

        Returns
        -------
        list
            Id list.
        """

        listId = np.asarray(listIdByQuery).squeeze().tolist()
        try:
            len(listId)
            return listId
        except TypeError:
            return [listId]

    def insertData(self, filterType, neighborStarMap):
        """Insert new star data into the local database.

        Parameters
        ----------
        filterType : FilterType
            Filter type.
        neighborStarMap : NbrStar
            Information of neighboring stars.
        """

        # List of bright star
        brightStarList = list(neighborStarMap.getId())

        # Check the existed bright star data based on ra and decl
        existIdList = []
        for ii in range(len(brightStarList)):
            raDec = neighborStarMap.getRaDecl()[brightStarList[ii]]
            if (self.searchRaDecl(filterType, raDec[0], raDec[1])):
                existIdList.append(brightStarList[ii])

        # Collect the lists not in database yet.
        # remainIDList is the bright star list. And allStarList is the list
        # contains the bright stars and related neighboring stars.
        remainIdList = []
        allStarList = []
        for id in brightStarList:
            if id not in existIdList:
                remainIdList.append(id)
                allStarList.append(id)
                for starID in neighborStarMap.getId()[id]:
                    # Make sure the starID is not in the allStarList yet
                    if starID not in allStarList:
                        allStarList.append(starID)

        # Insert the star data to local data base
        tableName = self._getTableName(filterType)
        for simobjID in allStarList:

            # Insert data
            command = "INSERT INTO " + tableName + " (simobjid, ra, decl, " + \
                      filterType.name.lower() + "mag, bright_star) " + \
                      "VALUES (?, ?, ?, ?, ?)"

            raDec = neighborStarMap.getRaDecl()[simobjID]
            mag = neighborStarMap.getMag(filterType)[simobjID]

            if simobjID in remainIdList:
                brightStar = True
            else:
                brightStar = False

            task = (int(simobjID), raDec[0], raDec[1], mag, brightStar)

            self.cursor.execute(command, task)

        # Commit the change to database
        self.connection.commit()

    def updateData(self, filterType, listID, listOfItemToChange,
                   listOfNewValue):
        """Update data based on the Id.

        Parameters
        ----------
        filterType : FilterType
            Filter type.
        listID : list[int]
            ID list to change.
        listOfItemToChange : list[str]
            Item list (simobjid, ra, decl, mag, bright_star) to change.
        listOfNewValue : list[int, float, or bool]
            New value list.

        Raises
        ------
        ValueError
            Not the correct type to update.
        """

        # Check the item can be updated or not
        for item in listOfItemToChange:
            if item not in ("simobjid", "ra", "decl", "mag", "bright_star"):
                raise ValueError("'%s' can not be updated." % item)

        # Update data based on the id
        tableName = self._getTableName(filterType)
        for ii in range(len(listID)):

            # Check the item is "mag" or not. If it is "mag", give the
            # related filter information.
            itemToChange = listOfItemToChange[ii]
            if (itemToChange == "mag"):
                itemToChange = filterType.name.lower() + itemToChange

            # Give the SQL command
            command = "UPDATE " + tableName + " SET " + \
                      itemToChange + "=" + str(listOfNewValue[ii]) + \
                      " WHERE id=?"
            self.cursor.execute(command, (listID[ii],))

        # Commit the change to database
        self.connection.commit()

    def deleteData(self, filterType, listID):
        """Delete data based on the Id.

        Parameters
        ----------
        filterType : FilterType
            Filter type
        listID : list[int]
            ID list to delete.
        """

        # Delete the data
        tableName = self._getTableName(filterType)
        for id in listID:
            command = "DELETE FROM " + tableName + " WHERE id=?"
            self.cursor.execute(command, (id,))

        # Commit the change to database
        self.connection.commit()

    def getAllId(self, filterType):
        """Get all ID in the database.

        Parameters
        ----------
        filterType : FilterType
            Filter type.

        Returns
        -------
        list
            ID list
        """

        # Print the table
        tableName = self._getTableName(filterType)
        command = "SELECT id FROM " + tableName
        self.cursor.execute(command)
        listIdByQuery = self.cursor.fetchall()

        return self._changeQueriedIdToList(listIdByQuery)


if __name__ == "__main__":
    pass
