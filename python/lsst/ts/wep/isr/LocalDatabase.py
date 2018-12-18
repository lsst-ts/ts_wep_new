import os, sqlite3
import numpy as np


class LocalDatabase(object):

    def __init__(self):

        self.connection = None
        self.cursor = None

    def connect(self, dbAdress):
        """
        
        Connect database based on the local path.
        
        Arguments:
            dbAdress {[string]} -- Path of local sqlite3 database.
        """

        self.connection = sqlite3.connect(dbAdress)
        self.cursor = self.connection.cursor()

    def printAll(self, tableName):
        """
        
        Print all data in table. This is only for the debug.
        
        Arguments:
            tableName {[string]} -- Name of table. For data butler, the tables are:
                                    "raw", "raw_visit", and "raw_skyTile" in 
                                    registry.sqlite3.
        """

        # Query all data in table
        command = "SELECT * FROM " + tableName
        self.cursor.execute(command)
        result = self.cursor.fetchall()

        # Print the table
        for line in result:
            print(line)

    def deleteData(self, visit, snap, raft, sensor, aFilter):
        """
        
        Delete data in tables: raw, raw_visit, and raw_skyTile.
        
        Arguments:
            visit {int} -- Visit time.
            snap {int} -- Snap time (0 or 1) means first/ second exposure.
            raft {[string]} -- Raft name.
            sensor {[string]} -- Sensor name.
            aFilter {[string]} -- Filter name (u, g, r, i, z, y).
        """

        # Query the "raw" database to get the dataID for each amplifier (channel)
        listId = self.__getRawId(visit, snap, raft, sensor, aFilter)

        # Delete the data in "raw" and "raw_skyTile"
        for dataId in listId:       
            # Delete the data in "raw"
            command = "DELETE FROM raw WHERE id=?"
            self.cursor.execute(command, (dataId,))

            # Delete the data in "raw_skyTile"
            command = "DELETE FROM raw_skyTile WHERE id=?"
            self.cursor.execute(command, (dataId,))

        # Delete the data in "raw_visit"
        command = "DELETE FROM raw_visit WHERE visit=?"
        self.cursor.execute(command, (visit,))            

        # Commit the change to database
        self.connection.commit()

    def insertData(self, visit, snap, raft, sensor, aFilter, taiObs, skyTile, expTime):
        """
        
        Insert data into tables: raw, raw_visit, and raw_skyTile.
        
        Arguments:
            visit {int} -- Visit time.
            snap {int} -- Snap time (0 or 1) means first/ second exposure.
            raft {[string]} -- Raft name.
            sensor {[string]} -- Sensor name.
            aFilter {[string]} -- Filter name (u, g, r, i, z, y).
            taiObs {[string]} -- International atomic time (TAI) of observation.
            skyTile {[int]} -- Sky tile.
            expTime {[float]} -- Exposure time.
        """

        # Amplifier list
        ampList = ["0", "0,0", "0,1", "0,2", "0,3", "0,4", "0,5", "0,6", "0,7", 
                   "1,0", "1,1", "1,2", "1,3", "1,4", "1,5", "1,6", "1,7"]

        # Update the condition for the WFS
        if (raft=="0,0" and sensor=="2,2") or (raft=="0,4" and sensor=="2,0") or (
            raft=="4,0" and sensor=="0,2") or (raft=="4,4" and sensor=="0,0"):
            ampList.append("1")

        # Check the data exists in "raw" table or not
        item = self.query("raw", visit=visit, snap=snap, raft=raft, sensor=sensor, aFilter=aFilter)

        if (len(item) == 0):
            # Insert new data into "raw" table
            command = "INSERT INTO raw (visit, filter, snap, raft, sensor, channel, taiObs, expTime) " + \
                      "VALUES (?, ?, ?, ?, ?, ?, ?, ?)"
            for channel in ampList:
                task = (visit, aFilter, snap, raft, sensor, channel, taiObs, expTime)
                self.cursor.execute(command, task)
                
            # Insert new data into "raw_visit" table
            command = "INSERT INTO raw_visit (visit, filter, taiObs, expTime) " + \
                      "VALUES (?, ?, ?, ?)"
            task = (visit, aFilter, taiObs, expTime)

            try:
                self.cursor.execute(command, task)
            except Exception as IntegrityError:
                print("Visit number %d exists in 'raw_visit' already." % visit)
            # self.cursor.execute(command, task)

            # Insert new data into "raw_skyTile" table
            # Query to get the primary key in "raw" table
            listId = self.__getRawId(visit, snap, raft, sensor, aFilter)
            for dataId in listId:
                command = "INSERT INTO raw_skyTile (id, skyTile) VALUES (?, ?)"
                task = (dataId, skyTile)
                self.cursor.execute(command, task)

            # Commit the change to database
            self.connection.commit()

            print("Add data into database.")

        else:
            print("Data exists in 'raw' table already.")

    def __getRawId(self, visit, snap, raft, sensor, aFilter):
        """
        
        Get the unique ID in "raw" table.
        
        Arguments:
            visit {int} -- Visit time (default: {0}).
            snap {int} -- Snap time (0 or 1) means first/ second exposure (default: {0}).
            raft {[string]} -- Raft name (default: {None}).
            sensor {[string]} -- Sensor name (default: {None}).
            aFilter {[string]} -- Filter name (u, g, r, i, z, y) (default: {None}).
        
        Returns:
            [int] -- Unique ID.
        """

        listId = []
        items = self.query("raw", visit=visit, snap=snap, raft=raft, sensor=sensor, aFilter=aFilter)
        for item in items:
            listId.append(item[0])

        return listId

    def disconnect(self):
        """
        Disconnect from the database.
        """

        self.cursor.close()
        self.connection.close()

    def query(self, tableName, dataId=0, visit=0, snap=0, raft=None, sensor=None, aFilter=None):
        """
        
        Query the database.
        
        Arguments:
            tableName {[string]} -- Name of table. For data butler, the tables are:
                                    "raw", "raw_visit", and "raw_skyTile" 
                                    in registry.sqlite3.
        
        Keyword Arguments:
            dataId {int} -- Unique data ID (default: {0}).
            visit {int} -- Visit time (default: {0}).
            snap {int} -- Snap time (0 or 1) means first/ second exposure (default: {0}).
            raft {[string]} -- Raft name (default: {None}).
            sensor {[string]} -- Sensor name (default: {None}).
            aFilter {[string]} -- Filter name (u, g, r, i, z, y) (default: {None}).
        
        Returns:
            [list] -- Queried data.
        """

        # Query the database of visit. 
        # Decide the command based on the tableName
        if (tableName == "raw"):
            command = "SELECT * FROM raw WHERE visit=%d AND snap=%d AND raft='%s' " + \
                      "AND sensor='%s' AND filter='%s'"
            query = command % (visit, snap, raft, sensor, aFilter)
        elif (tableName == "raw_skyTile"):
            command = "SELECT * FROM raw_skyTile WHERE id='%s'"
            query = command % dataId
        elif (tableName == "raw_visit"):
            command = "SELECT * FROM raw_visit WHERE visit=%d and filter='%s'"
            query = command % (visit, aFilter)
        else:
            print("No such database. Only 'raw', 'raw_visit', and 'raw_skyTile'.")

        # Execute the query
        self.cursor.execute(query)

        return self.cursor.fetchall()


if __name__ == "__main__":
    pass
