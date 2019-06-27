from lsst.ts.wep.Utility import CentroidFindType
from lsst.ts.wep.cwfs.CentroidRandomWalk import CentroidRandomWalk
from lsst.ts.wep.cwfs.CentroidOtsu import CentroidOtsu


class CentroidFindFactory(object):
    """Factory for creating the centroid find object to calculate the centroid
    of donut."""

    @staticmethod
    def createCentroidFind(centroidFindType):
        """Create the centroid find object.

        Parameters
        ----------
        centroidFindType : CentroidFindType
            Algorithm to find the centroid.

        Returns
        -------
        CentroidRandomWalk, CentroidOtsu
            Centroid find object.

        Raises
        ------
        ValueError
            The centroid find type is not supported.
        """

        if (centroidFindType == CentroidFindType.RandomWalk):
            return CentroidRandomWalk()
        elif (centroidFindType == CentroidFindType.Otsu):
            return CentroidOtsu()
        else:
            raise ValueError("The %s is not supported." % centroidFindType)


if __name__ == "__main__":
    pass
