# This file is part of ts_wep.
#
# Developed for the LSST Telescope and Site Systems.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np

from lsst.ts.wep.bsc.LocalDatabase import LocalDatabase


class LocalDatabaseForStarFile(LocalDatabase):

    PRE_TABLE_NAME = "StarTable"

    def createTable(self, filterType):
        """Create the table in database.

        Parameters
        ----------
        filterType : FilterType
            Filter type.

        Raises
        ------
        ValueError
            Table exists in database already.
        """

        tableName = self._getTableName(filterType)
        if self._tableIsInDb(tableName):
            raise ValueError("%s exists in database already." % tableName)

        # Create the table
        command = (
            "CREATE TABLE %s "
            "(id INTEGER PRIMARY KEY, simobjid INTEGER NOT NULL, "
            "ra REAL, decl REAL, %smag REAL, bright_star NUMERIC)"
        ) % (tableName, filterType.name.lower())

        self.cursor.execute(command)

        # Commit the change to database
        self.connection.commit()

    def _tableIsInDb(self, tableName):
        """Check the specific table exists in the database or not.

        Parameters
        ----------
        tableName : str
            Table name.

        Returns
        -------
        bool
            Table exists or not.
        """

        # Query the specific table name
        command = "PRAGMA table_info(%s)" % tableName
        self.cursor.execute(command)
        data = self.cursor.fetchall()

        # Check the table exists or not
        if len(data) == 0:
            return False
        else:
            return True

    def insertDataByFile(self, skyFilePath, filterType, skiprows=1):
        """Insert the sky data by file.

        Parameters
        ----------
        skyFilePath : str
            Sky data file path.
        filterType : FilterType
            Filter type.
        skiprows : int, optional
            Skip the first 'skiprows' lines. (the default is 1.)
        """

        # Get the data
        skyData = np.loadtxt(skyFilePath, skiprows=skiprows)

        # Only consider the non-empty data
        if len(skyData) != 0:

            # Change to 2D array if the input is 1D array
            if skyData.ndim == 1:
                skyData = np.expand_dims(skyData, axis=0)

            # Add the star
            tableName = self._getTableName(filterType)
            for ii in range(len(skyData)):
                # Insert data
                command = (
                    "INSERT INTO %s "
                    "(simobjid, ra, decl, %smag, bright_star) "
                    "VALUES (?, ?, ?, ?, ?)"
                ) % (tableName, filterType.name.lower())

                simobjID, ra, decl, mag = skyData[ii]
                task = (int(simobjID), ra, decl, mag, 0)

                self.cursor.execute(command, task)

            # Commit the change to database
            self.connection.commit()

    def deleteTable(self, filterType):
        """Delete the table in database.

        Parameters
        ----------
        filterType : FilterType
            Filter type.
        """

        # Delete the table
        tableName = self._getTableName(filterType)
        command = "DROP TABLE IF EXISTS %s" % tableName
        self.cursor.execute(command)

        # Commit the change to database
        self.connection.commit()


if __name__ == "__main__":
    pass
