from lsst.ts.wep.bsc.LocalDatabase import LocalDatabase
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
        LocalDatabase
            BSC database object.

        Raises
        ------
        ValueError
            The database type does not match.
        """

        if (bscDbType == BscDbType.LocalDb):
            return LocalDatabase()
        else:
            raise ValueError("The database type does not match.")


if __name__ == "__main__":
    pass
