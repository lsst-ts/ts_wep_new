import numpy as np


class NeighboringStar(object):

    def __init__(self):
        """Initialize the neighboring star class."""

        self.SimobjID = {} 
        self.RaDecl = {}
        self.RaDeclInPixel = {}

        self.LSSTMagU = {}
        self.LSSTMagG = {}
        self.LSSTMagR = {}
        self.LSSTMagI = {}
        self.LSSTMagZ = {}
        self.LSSTMagY = {}

    def addStar(self, stars, indexCandidate, indexNeighboringStar,
                cameraFilter):
        """Add the information of neightboring stars based on the candidate
        star.

        Parameters
        ----------
        stars : StarData
            Star information.
        indexCandidate : int
            Index of candidate stars.
        indexNeighboringStar : int
            Index of neighboring stars of specific candidate star.
        cameraFilter : str
            Filter type of camera: u, g, r, i, z, y.
        """

        self.SimobjID[stars.SimobjID[indexCandidate]] = np.array(
                            stars.SimobjID)[[indexNeighboringStar]].tolist()

        # Collect coordinates and magnitude of candidate and neighboring stars 
        indexStar = np.append(indexNeighboringStar, indexCandidate)
        for index in indexStar:
            self.RaDecl[stars.SimobjID[index]] = (stars.RA[index],
                                                  stars.Decl[index])
            self.RaDeclInPixel[stars.SimobjID[index]] = \
                        (stars.RAInPixel[index], stars.DeclInPixel[index])

            if (cameraFilter == stars.FilterU):
                self.LSSTMagU[stars.SimobjID[index]] = stars.LSSTMagU[index]

            elif (cameraFilter == stars.FilterG):
                self.LSSTMagG[stars.SimobjID[index]] = stars.LSSTMagG[index]

            elif (cameraFilter == stars.FilterR):
                self.LSSTMagR[stars.SimobjID[index]] = stars.LSSTMagR[index]

            elif (cameraFilter == stars.FilterI):
                self.LSSTMagI[stars.SimobjID[index]] = stars.LSSTMagI[index]

            elif (cameraFilter == stars.FilterZ):
                self.LSSTMagZ[stars.SimobjID[index]] = stars.LSSTMagZ[index]

            elif (cameraFilter == stars.FilterY):
                self.LSSTMagY[stars.SimobjID[index]] = stars.LSSTMagY[index]


if __name__ == "__main__":
    pass
