import numpy as np

from lsst.ts.wep.Utility import FilterType
from lsst.ts.wep.bsc.StarData import StarData


class DefaultDatabase(object):

    # Value to decide the query will cross the RA=0 or not
    STD_DEV_SPLIT = 20.0

    def __init__(self):
        """Initialize the default database class."""

        self.connection = None
        self.cursor = None

    def disconnect(self):
        """Disconnect the database."""

        self.cursor.close()
        self.connection.close()

    def query(self, filterType, corner1, corner2, corner3, corner4):
        """Query the database for stars within an area.

        Parameters
        ----------
        filterType : FilterType
            Filter type.
        corner1 : tuple
            The first corner of the sensor defined as (RA, Decl).
        corner2 : tuple
            The second corner of the sensor defined as (RA, Decl).
        corner3 : tuple
            The third corner of the sensor defined as (RA, Decl).
        corner4 : tuple
            The fourth corner of the sensor defined as (RA, Decl).

        Returns
        ----------
        StarData
            Star information.
        """

        ra = [corner1[0], corner2[0], corner3[0], corner4[0]]
        decl = [corner1[1], corner2[1], corner3[1], corner4[1]]
        top = max(decl)
        bottom = min(decl)
        left = min(ra)
        right = max(ra)

        # Need to change this query method that divides the area with 2 parts.
        raStddev = np.std(ra)

        # Query regions crosses the RA = 0
        if (raStddev >= self.STD_DEV_SPLIT):

            # Query the left and right regions
            left = max([x for x in ra if x < 180])
            right = min([x for x in ra if x >= 180])
            above0Set = self._queryTable(filterType, top, bottom, 0, left)
            below0Set = self._queryTable(filterType, top, bottom, right, 360)

            # Combine the query results
            starId = np.append(above0Set.getId(), below0Set.getId())
            starRa = np.append(above0Set.getRA(), below0Set.getRA())
            starDecl = np.append(above0Set.getDecl(), below0Set.getDecl())

            lsstMagU = np.append(above0Set.getMag(FilterType.U),
                                 below0Set.getMag(FilterType.U))
            lsstMagG = np.append(above0Set.getMag(FilterType.G),
                                 below0Set.getMag(FilterType.G))
            lsstMagR = np.append(above0Set.getMag(FilterType.R),
                                 below0Set.getMag(FilterType.R))
            lsstMagI = np.append(above0Set.getMag(FilterType.I),
                                 below0Set.getMag(FilterType.I))
            lsstMagZ = np.append(above0Set.getMag(FilterType.Z),
                                 below0Set.getMag(FilterType.Z))
            lsstMagY = np.append(above0Set.getMag(FilterType.Y),
                                 below0Set.getMag(FilterType.Y))

            return StarData(starId, starRa, starDecl, lsstMagU, lsstMagG,
                            lsstMagR, lsstMagI, lsstMagZ, lsstMagY)

        # Query regions does not cross the RA = 0
        else:
            return self._queryTable(filterType, top, bottom, left, right)

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

        Raises
        ------
        NotImplementedError
            Child class should implemented this.
        """

        raise NotImplementedError("Child class should implemented this.")

    def connect(self):
        """Connect the database.

        The child class needs to concrete the connection and cursor attributes
        with the Connection and Cursor objects. For example, the SQLite
        Connection and Cursor objects.

        Raises
        ------
        NotImplementedError
            Child class should implemented this.
        """

        raise NotImplementedError("Child class should implemented this.")


if __name__ == "__main__":
    pass
