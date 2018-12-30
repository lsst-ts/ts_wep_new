from lsst.ts.wep.bsc.LocalDatabase import LocalDatabase
from lsst.ts.wep.Utility import BscDbType


class DatabaseFactory(object):

    @staticmethod
    def createDb(dbType):
        """Create the database.

        Parameters
        ----------
        dbType : BscDbType
            Database type.

        Returns
        -------
        LocalDatabase
            Database object.

        Raises
        ------
        ValueError
            The database type does not match.
        """

        if (dbType == BscDbType.LocalDb):
            return LocalDatabase()
        else:
            raise ValueError("The database type does not match.")


if __name__ == "__main__":
    pass
