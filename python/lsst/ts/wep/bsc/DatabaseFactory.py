from lsst.ts.wep.bsc.LocalDatabase import LocalDatabase
from lsst.ts.wep.bsc.LocalDatabaseForStarFile import LocalDatabaseForStarFile

from lsst.ts.wep.Utility import BscDbType


class DatabaseFactory(object):

    @staticmethod
    def createDb(bscDbType):
        """Create the database.

        Parameters
        ----------
        bscDbType : BscDbType
            Bright star catalog (BSC) database type.

        Returns
        -------
        LocalDatabase or LocalDatabaseForStarFile
            BSC database object.

        Raises
        ------
        ValueError
            The database type does not match.
        """

        if (bscDbType == BscDbType.LocalDb):
            return LocalDatabase()
        elif (bscDbType == BscDbType.LocalDbForStarFile):
            return LocalDatabaseForStarFile()
        else:
            raise ValueError("The database type does not match.")


if __name__ == "__main__":
    pass
