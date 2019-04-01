import numpy as np
from scipy.spatial.distance import cdist

from lsst.ts.wep.Utility import FilterType
from lsst.ts.wep.bsc.NbrStar import NbrStar


class StarData(object):

    def __init__(self, starId, ra, decl, lsstMagU, lsstMagG, lsstMagR,
                 lsstMagI, lsstMagZ, lsstMagY):
        """Initialize the star data class.

        Parameters
        ----------
        starId : int, list[int], or numpy.ndarray[int]
            Star Id.
        ra : float, list[float], or numpy.ndarray[float]
            Star right ascension in degree.
        decl : float, list[float], or numpy.ndarray[float]
            Star declination in degree.
        lsstMagU : float, list[float], or numpy.ndarray[float]
            Star magnitude under the LSST u filter band.
        lsstMagG : float, list[float], or numpy.ndarray[float]
            Star magnitude under the LSST g filter band.
        lsstMagR : float, list[float], or numpy.ndarray[float]
            Star magnitude under the LSST r filter band.
        lsstMagI : float, list[float], or numpy.ndarray[float]
            Star magnitude under the LSST i filter band.
        lsstMagZ : float, list[float], or numpy.ndarray[float]
            Star magnitude under the LSST z filter band.
        lsstMagY : float, list[float], or numpy.ndarray[float]
            Star magnitude under the LSST y filter band.
        """

        self.detector = ""
        self.starId = np.array([], dtype=int)
        self.ra = np.array([])
        self.decl = np.array([])

        self.raInPixel = np.array([])
        self.declInPixel = np.array([])

        self.lsstMagU = np.array([])
        self.lsstMagG = np.array([])
        self.lsstMagR = np.array([])
        self.lsstMagI = np.array([])
        self.lsstMagZ = np.array([])
        self.lsstMagY = np.array([])

        self.setId(starId)
        self.setRA(ra)
        self.setDecl(decl)

        self.setMag(FilterType.U, lsstMagU)
        self.setMag(FilterType.G, lsstMagG)
        self.setMag(FilterType.R, lsstMagR)
        self.setMag(FilterType.I, lsstMagI)
        self.setMag(FilterType.Z, lsstMagZ)
        self.setMag(FilterType.Y, lsstMagY)

    def setId(self, starId):
        """Set the star Id.

        Parameters
        ----------
        starId : int, list[int], or numpy.ndarray[int]
            Star Id.
        """

        starIdArray = self._changeToNpArrayIfNeeded(starId)
        self.starId = starIdArray.astype(int)

    def setRA(self, ra):
        """Set the star right ascension (RA) in degree.

        Parameters
        ----------
        ra : float, list[float], or numpy.ndarray[float]
            Star right ascension in degree.
        """

        self.ra = self._changeToNpArrayIfNeeded(ra)

    def setDecl(self, decl):
        """Set the star declination (Decl) in degree.

        Parameters
        ----------
        decl : float, list[float], or numpy.ndarray[float]
            Star declination in degree.
        """

        self.decl = self._changeToNpArrayIfNeeded(decl)

    def setMag(self, filterType, mag):
        """Set the star magnitude.

        Parameters
        ----------
        filterType : FilterType
            Filter type.
        mag : float, list[float], or numpy.ndarray[float]
            Star magnitude under the specific filter type.

        Raises
        ------
        ValueError
            No filter type matches.
        """

        if (filterType == FilterType.U):
            self.lsstMagU = self._changeToNpArrayIfNeeded(mag)
        elif (filterType == FilterType.G):
            self.lsstMagG = self._changeToNpArrayIfNeeded(mag)
        elif (filterType == FilterType.R):
            self.lsstMagR = self._changeToNpArrayIfNeeded(mag)
        elif (filterType == FilterType.I):
            self.lsstMagI = self._changeToNpArrayIfNeeded(mag)
        elif (filterType == FilterType.Z):
            self.lsstMagZ = self._changeToNpArrayIfNeeded(mag)
        elif (filterType == FilterType.Y):
            self.lsstMagY = self._changeToNpArrayIfNeeded(mag)
        else:
            raise ValueError("No filter type matches.")

    def _changeToNpArrayIfNeeded(self, val):
        """Change the value type to the numpy array if it is needed.

        Parameters
        ----------
        val : int, float, list, or numpy.ndarray
            Value.

        Returns
        -------
        numpy.ndarray
            Value in the numpy array data type.
        """

        if isinstance(val, np.ndarray):
            return val
        elif isinstance(val, list):
            return np.array(val)
        else:
            return np.array([val])

    def getId(self):
        """Get the star Id.

        Returns
        -------
        numpy.ndarray
            Star Id.
        """

        return self.starId

    def getRA(self):
        """Get the star right ascension (RA) in degree.

        Returns
        -------
        numpy.ndarray
            Star RA in degree.
        """

        return self.ra

    def getDecl(self):
        """Get the star declination in degree.

        Returns
        -------
        numpy.ndarray
            Star declination in degree.
        """

        return self.decl

    def getRaInPixel(self):
        """Get the star right ascension (RA) in pixel.

        Returns
        -------
        numpy.ndarray
            Star RA in pixel.
        """

        return self.raInPixel

    def getDeclInPixel(self):
        """Get the star declination in pixel.

        Returns
        -------
        numpy.ndarray
            Star declination in pixel.
        """

        return self.declInPixel

    def getMag(self, filterType):
        """Get the star magnitude.

        Parameters
        ----------
        filterType : FilterType
            Filter type.

        Returns
        -------
        numpy.ndarray
            Star magnitude.

        Raises
        ------
        ValueError
            No filter type matches.
        """

        if (filterType == FilterType.U):
            return self.lsstMagU
        elif (filterType == FilterType.G):
            return self.lsstMagG
        elif (filterType == FilterType.R):
            return self.lsstMagR
        elif (filterType == FilterType.I):
            return self.lsstMagI
        elif (filterType == FilterType.Z):
            return self.lsstMagZ
        elif (filterType == FilterType.Y):
            return self.lsstMagY
        else:
            raise ValueError("No filter type matches.")

    def getDetector(self):
        """Get the detector.

        Returns
        -------
        str
            Detector.
        """

        return self.detector

    def setDetector(self, detector):
        """Set the detector name for this set of stars.

        Parameters
        ----------
        detector : str
            The name of the detector.
        """

        self.detector = detector

    def setRaInPixel(self, raInPixel):
        """Set the RA pixel data for this set of stars.

        Parameters
        ----------
        raInPixel : float, list[float], or numpy.array[float]
            The ra pixel coordinate each star falls on the detector.
        """

        self.raInPixel = self._changeToNpArrayIfNeeded(raInPixel)

    def setDeclInPixel(self, declInPixel):
        """Set the Decl pixel data for this set of stars.

        Parameters
        ----------
        declInPixel : float, list[float], or numpy.array[float]
            The decl pixel coordinate each star falls on the detector.
        """

        self.declInPixel = self._changeToNpArrayIfNeeded(declInPixel)

    def checkCandidateStars(self, filterType, lowMag, highMag):
        """Check the candidate stars based on the magnitude.

        Parameters
        ----------
        filterType : FilterType
            Filter type.
        lowMag : float
            Lower boundary of magnitude.
        highMag : float
            Higher boundary of magnitude.

        Returns
        -------
        list[int]
            List of index candidate.
        """

        if (len(self.ra) != 0):

            if (filterType == FilterType.U):
                valArray = self.lsstMagU

            elif (filterType == FilterType.G):
                valArray = self.lsstMagG

            elif (filterType == FilterType.R):
                valArray = self.lsstMagR

            elif (filterType == FilterType.I):
                valArray = self.lsstMagI

            elif (filterType == FilterType.Z):
                valArray = self.lsstMagZ

            elif (filterType == FilterType.Y):
                valArray = self.lsstMagY

            idxCand = self._getIdxCandByBndry(valArray, lowMag, highMag)

        else:
            idxCand = []

        return idxCand

    def _getIdxCandByBndry(self, valArray, lowMag, highMag):
        """Get the index candidate list based on the boundary of magnitude.

        Parameters
        ----------
        valArray : list or numpy.ndarray
            Value array
        lowMag : float
            Lower boundary of magnitude.
        highMag : float
            Higher boundary of magnitude.

        Returns
        -------
        list
            List of index candidate.
        """

        idxCand = [idx for idx in range(len(valArray))
                   if lowMag <= valArray[idx] <= highMag]

        return idxCand

    def getNeighboringStar(self, idxCand, maxDist, filterType,
                           maxNumOfNbrStar=0):
        """Get the neighboring stars of candidate stars based on the specific
        max distance and number of neighboring star.

        Parameters
        ----------
        idxCand : list[int]
            List of index of candidate star in "stars" data.
        maxDist : float
            Maximum distance in pixel.
        filterType : FilterType
            Filter type.
        maxNumOfNbrStar : int, optional
            Maximum number of neighboring star for each candidate star (the
            default is 0)

        Returns
        -------
        NbrStar
            Information of neighboring stars.
        """

        nbrStar = NbrStar()

        # Calculate the distance in pixel between candidate stars and all stars
        numOfIdxCand = len(idxCand)
        if (numOfIdxCand != 0):
            allStarXY = np.array([self.raInPixel,
                                  self.declInPixel]).transpose()
            candidateStarXY = allStarXY[np.array(idxCand), :]
            starDistances = cdist(candidateStarXY, allStarXY)

            # Decide the number of neighboring stars
            for ii in range(numOfIdxCand):
                idxNbrStar = np.where(starDistances[ii, :] < maxDist)[0]

                # Delete candidate star itself
                idxNbrStar = np.delete(idxNbrStar,
                                       np.where(idxNbrStar == idxCand[ii]))

                # Remove the candidate star if there is the neighboring star
                # brighter than itself
                if (filterType == FilterType.U):
                    magSelf = self.lsstMagU[idxCand[ii]]
                    magNbrStar = self.lsstMagU[[idxNbrStar]]

                elif (filterType == FilterType.G):
                    magSelf = self.lsstMagG[idxCand[ii]]
                    magNbrStar = self.lsstMagG[[idxNbrStar]]

                elif (filterType == FilterType.R):
                    magSelf = self.lsstMagR[idxCand[ii]]
                    magNbrStar = self.lsstMagR[[idxNbrStar]]

                elif (filterType == FilterType.I):
                    magSelf = self.lsstMagI[idxCand[ii]]
                    magNbrStar = self.lsstMagI[[idxNbrStar]]

                elif (filterType == FilterType.Z):
                    magSelf = self.lsstMagZ[idxCand[ii]]
                    magNbrStar = self.lsstMagZ[[idxNbrStar]]

                elif (filterType == FilterType.Y):
                    magSelf = self.lsstMagY[idxCand[ii]]
                    magNbrStar = self.lsstMagY[[idxNbrStar]]

                if ((np.where(magNbrStar < magSelf)[0]).size != 0):
                    brighterNeighbor = True
                else:
                    brighterNeighbor = False

                # Restrict the maximum number of neighboring stars
                if (len(idxNbrStar) > maxNumOfNbrStar):
                    highNeighboringStar = True
                else:
                    highNeighboringStar = False

                # Record the information of neighboring stars
                if (brighterNeighbor is False and highNeighboringStar is False):
                    nbrStar.addStar(self, idxCand[ii], idxNbrStar, filterType)

        return nbrStar


if __name__ == "__main__":
    pass
