import os, unittest
import numpy as np

from bsc.LocalDatabase import LocalDatabase

from lsst.ts.wep.Utility import getModulePath

class LocalDatabaseDecorator(LocalDatabase):

    def createTable(self, aFilter, tableName):
        """
        
        Create the table in database.
        
        Arguments:
            aFilter {[str]} -- A filter type ("u", "g", "r", "i", "z", "y").
            tableName {[str]} -- Table name.
        """

        # Check the filter type
        self.__checkFilterType(aFilter)

        # Check the table is in the database or not
        if (self.checkTableInDb(tableName)):
            raise RuntimeError("The table:'%s' exists in database already." % tableName)
        
        # Create the table
        command = "CREATE TABLE %s" % tableName
        command += "(id INTEGER PRIMARY KEY, simobjid INTEGER NOT NULL, "
        command += "ra REAL, decl REAL, %smag REAL, bright_star NUMERIC)" % aFilter
        self.cursor.execute(command)

        # Commit the change to database
        self.connection.commit()

    def insertDataByFile(self, aFilter, tableName, skyFilePath, skiprows=1):
        """
        
        Insert the sky data by file
        
        Arguments:
            aFilter {[str]} -- A filter type ("u", "g", "r", "i", "z", "y").
            tableName {[str]} -- Table name.
            skyFilePath {[str]} -- Sky data file path.
        
        Keyword Arguments:
            skiprows {int} -- Skip the first "skiprows" lines. (default: {1})
        """

        # Check the filter type
        self.__checkFilterType(aFilter)

        # Get the data
        skyData = np.loadtxt(skyFilePath, skiprows=skiprows)

        # Add the star
        command = "Insert"
        for ii in range(len(skyData)):

            # Insert data
            command = "INSERT INTO " + tableName + \
                      " (simobjid, ra, decl, " + aFilter + "mag, bright_star) " + \
                      "VALUES (?, ?, ?, ?, ?)"

            simobjID, ra, decl, mag = skyData[ii]

            task = (int(simobjID), ra, decl, mag, 0)

            self.cursor.execute(command, task)

        # Commit the change to database
        self.connection.commit()

    def deleteTable(self, tableName):
        """
        
        Delete the table in database.
        
        Arguments:
            tableName {[str]} -- Table name.
        """
        
        # Delete the table
        command = "DROP TABLE IF EXISTS %s" % tableName
        self.cursor.execute(command)

        # Commit the change to database
        self.connection.commit()

    def __checkFilterType(self, aFilter):
        """
        
        Check the filter type.
        
        Arguments:
             aFilter {[str]} -- A filter type ("u", "g", "r", "i", "z", "y").
        
        Raises:
            ValueError -- Filter does not match.
        """

        if aFilter not in (self.FilterU, self.FilterG, self.FilterR, self.FilterI, 
                            self.FilterZ, self.FilterY):
            raise ValueError("The filter type: '%s' is not allowed." % aFilter)

    def checkTableInDb(self, tableName):
        """
        
        Check the specific table exists in the database or not.
        
        Arguments:
            tableName {[str]} -- Table name.
        
        Returns:
            [bool] -- Table exists or not.
        """
        
        # Query the specific table name
        command = "PRAGMA table_info(%s)" % tableName
        self.cursor.execute(command)
        data = self.cursor.fetchall()

        # Check the table exists or not
        if (len(data) == 0):
            isExist = False
        else:
            isExist = True

        return isExist

class LocalDatabaseDecoratorTest(unittest.TestCase):
    """
    Test the LocalDatabaseDecorator. 
    """

    def setUp(self):

        # Get the path of module
        self.modulePath = getModulePath()

        # Address of local database
        dbAdress = os.path.join(self.modulePath, "test", "bsc.db3")

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
        skyFilePath = os.path.join(self.modulePath, "test", "skyComCamInfo.txt")
        self.db.insertDataByFile(aFilter, tableName, skyFilePath)

        self.db.deleteTable(tableName)
        self.assertFalse(self.db.checkTableInDb(tableName))

if __name__ == "__main__":

    # Do the unit test
    unittest.main()