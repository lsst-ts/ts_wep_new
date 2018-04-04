import unittest

class DefocalImage(object):

    def __init__(self, intraImg=None, extraImg=None):
        """
        
        Initialize the DefocalImage class.
        
        Keyword Arguments:
            intraImg {[ndarray]} -- Intra-defocal image. (default: {None})
            extraImg {[ndarray]} -- Extra-defocal image. (default: {None})
        """
        
        # Defocal images
        self.intraImg = intraImg
        self.extraImg = extraImg

    def setImg(self, intraImg=None, extraImg=None):
        """
        
        Set the image.
        
        Keyword Arguments:
            intraImg {[ndarray]} -- Intra-defocal image. (default: {None})
            extraImg {[ndarray]} -- Extra-defocal image. (default: {None})
        """

        if (intraImg is not None):
            self.intraImg = intraImg

        if (extraImg is not None):
            self.extraImg = extraImg

class DonutImage(DefocalImage):

    def __init__(self, starId, pixelX, pixelY, fieldX, fieldY, intraImg=None, extraImg=None):
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

class DefocalImageTest(unittest.TestCase):

    """ 
    Test the function of DefocalImage.
    """

    def setUp(self):
        
        self.defocalImg = DefocalImage()

    def testFunction(self):
        
        intraImg = 1
        extraImg = 2
        self.defocalImg.setImg(intraImg=intraImg, extraImg=extraImg)
        self.assertEqual(self.defocalImg.intraImg, intraImg)
        self.assertEqual(self.defocalImg.extraImg, extraImg)

class DonutImageTest(unittest.TestCase):

    """ 
    Test the function of DonutImage.
    """

    def setUp(self):
        
        starId = 0 
        pixelX = 1 
        pixelY = 1 
        fieldX = 2 
        fieldY = 2
        self.donutImg = DonutImage(starId, pixelX, pixelY, fieldX, fieldY)

    def testFunction(self):
        
        self.assertEqual(self.donutImg.pixelX, self.donutImg.pixelY)

        zer4UpNm = 10
        self.donutImg.setWfErr(zer4UpNm)
        self.assertEqual(self.donutImg.zer4UpNm, zer4UpNm)

if __name__ == "__main__":
    
    # Do the unit test
    unittest.main()