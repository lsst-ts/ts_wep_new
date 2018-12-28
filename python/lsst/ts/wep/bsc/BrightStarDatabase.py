import pymssql
from numpy import std
from decimal import Decimal

from lsst.ts.wep.bsc.Filter import Filter
from lsst.ts.wep.bsc.StarData import StarData


class BrightStarDatabase(object):

    def __init__(self):
        self.connection = None
        self.cursor = None

        # Get the filter type defined in the Filter class
        afilter = Filter()
        self.FilterU = afilter.FilterU
        self.FilterG = afilter.FilterG
        self.FilterR = afilter.FilterR
        self.FilterI = afilter.FilterI
        self.FilterZ = afilter.FilterZ
        self.FilterY = afilter.FilterY

        # Value to decide the query will cross the RA=0 or not
        self.stddevSplit = 20.0
    
    def connect(self, host, user, password, database):
        """
        
        Connects to the database.
        
        Arguments:
            host {[string]} -- The host name / ip and port of the database.
            user {[string]} -- The user name to connect as.
            password {[string]} -- The password for the user.
            database {[string]} -- The database to use.
        """
        self.connection = pymssql.connect(host, user, password, database)
        self.cursor = self.connection.cursor()

    def query(self, tableName, cameraFilter, corner1, corner2, corner3, corner4):
        """
        
        Queries the database for stars within an area.
        
        Arguments:
            tableName {[string]} -- Table name in database.
            cameraFilter {[string]} -- Filter type of camera: u, g, r, i, z, y.
            corner1 {[float]} -- The first corner of the sensor defined as (RA, Decl).
            corner2 {[float]} -- The second corner of the sensor defined as (RA, Decl).
            corner3 {[float]} -- The third corner of the sensor defined as (RA, Decl).
            corner4 {[float]} -- The fourth corner of the sensor defined as (RA, Decl).
        
        Returns:
            [metadata] -- Do the database query.
        """

        ra = [corner1[0], corner2[0], corner3[0], corner4[0]]
        decl = [corner1[1], corner2[1], corner3[1], corner4[1]]
        top = max(decl)
        bottom = min(decl)
        left = min(ra)
        right = max(ra)
        
        # Need to change this query method that divides the area with 2 parts.
        # Use the stddevSplit might not be a good idea.
        raStddev = std(ra)
        if raStddev >= self.stddevSplit:
            left = max([x for x in ra if x < 180])
            right = min([x for x in ra if x >= 180])
            above0Set = self.queryTable(tableName, cameraFilter, top, bottom, 0, left)
            below0Set = self.queryTable(tableName, cameraFilter, top, bottom, right, 360)
            return StarData(above0Set.SimobjID + below0Set.SimobjID,
                            above0Set.RA + below0Set.RA, 
                            above0Set.Decl + below0Set.Decl, 
                            above0Set.LSSTMagU + below0Set.LSSTMagU, 
                            above0Set.LSSTMagG + below0Set.LSSTMagG, 
                            above0Set.LSSTMagR + below0Set.LSSTMagR, 
                            above0Set.LSSTMagI + below0Set.LSSTMagI, 
                            above0Set.LSSTMagZ + below0Set.LSSTMagZ, 
                            above0Set.LSSTMagY + below0Set.LSSTMagY)
        else:                            
            return self.queryTable(tableName, cameraFilter, top, bottom, left, right)

    def queryTable(self, tableName, cameraFilter, top, bottom, left, right):
        """
        
        Queries the database for stars within an area.
                
        Arguments:
            tableName {[string]} -- Table name in database.
            cameraFilter {[string]} -- Filter type of camera: u, g, r, i, z, y.
            top {[float]} -- The top edge of the box (Decl).
            bottom {[float]} -- The bottom edge of the box (Decl).
            left {[float]} -- The left edge of the box (RA).
            right {[float]} -- The right edge of the box (RA).
        
        Returns:
            [StarData] -- Star information.
        """


        # Get the lsst filter magnitudes of stars in a certain range
        command = "SELECT simobjid, ra, decl, " + cameraFilter + "mag " + \
                  "FROM " + tableName + \
                  " WHERE decl <= %f AND decl >= %f AND ra >= %f AND ra <= %f"
                
        query = command % (top, bottom, left, right)

        self.cursor.execute(query)

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
            
            # It is noted that the data type of simobjid is big interger in UW database 
            simobjid.append(item[0])
            ra.append(item[1])
            decl.append(item[2])
           
            if (cameraFilter == self.FilterU):              
                lsstMagU.append(item[3])    

            elif (cameraFilter == self.FilterG):              
                lsstMagG.append(item[3])

            elif (cameraFilter == self.FilterR):
                lsstMagR.append(item[3])

            elif (cameraFilter == self.FilterI):
                lsstMagI.append(item[3])

            elif (cameraFilter == self.FilterZ):
                lsstMagZ.append(item[3])

            elif (cameraFilter == self.FilterY):
                lsstMagY.append(item[3])

        return StarData(simobjid, ra, decl, lsstMagU, lsstMagG, lsstMagR,  
                        lsstMagI, lsstMagZ, lsstMagY)                

    def disconnect(self):
        """
        Disconnects from the database.
        """
        self.cursor.close()
        self.connection.close()

    def searchRaDecl(self, tableName, ra, decl):
        """
        
        Search the star simobjid based on ra, decl.
        
        Arguments:
            tableName {[string]} -- Table name in database.
            ra {[float]} -- ra in degree (0 deg - 360 deg).
            decl {[float]} -- decl in degree (-90 deg - 90 deg).
        
        Returns:
            [decimal] -- Star simobjid in remote UW databse.
        """

        # Compare ra and decl to see the existance of star in database
        command = "SELECT simobjid FROM " + tableName + " WHERE ra = %f AND decl = %f" 
        query = command % (ra, decl)
        self.cursor.execute(query)
        
        return self.cursor.fetchall()


if __name__ == "__main__":
    pass
