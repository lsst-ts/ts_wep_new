from lsst.ts.wep.DefocalImage import DefocalImage


class DonutImage(DefocalImage):

    def __init__(self, starId, pixelX, pixelY, fieldX, fieldY, intraImg=None,
                 extraImg=None):
        """
        
        Initialize the DonutImage class.
        
        Arguments:
            starId {[int]} -- star ID.
            pixelX {[float]} -- Pixel x.
            pixelY {[float]} -- Pixel y.
            fieldX {[float]} -- Field x in degree.
            fieldY {[float]} -- Field y in degree.
        
        Keyword Arguments:
            intraImg {[ndarray]} -- Intra-defocal image. (default: {None})
            extraImg {[ndarray]} -- Extra-defocal image. (default: {None})
        """

        super(DonutImage, self).__init__(intraImg=intraImg, extraImg=extraImg)
        self.starId = int(starId)
        self.pixelX = pixelX
        self.pixelY = pixelY
        self.fieldX = fieldX
        self.fieldY = fieldY

        # Wavefront eror in annular Zk in nm (z4-z22)
        self.zer4UpNm = None

    def setWfErr(self, zer4UpNm):
        """
        
        Set the wavefront error in annular Zk in nm (z4-z22).
        
        Arguments:
            zer4UpNm {[ndarray]} -- z4 to z22 in nm.
        """

        self.zer4UpNm = zer4UpNm


if __name__ == "__main__":
    pass
