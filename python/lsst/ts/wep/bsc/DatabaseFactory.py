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
        """

        if (dbType == BscDbType.LocalDb):
            return LocalDatabase()


if __name__ == "__main__":
    pass
