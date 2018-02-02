class DefocalImage(object):

    def __init__(self, intraImg=None, extraImg=None):
        """
        
        Initialize the DefocalImage class.
        
        Keyword Arguments:
            intraImg {[ndarray]} -- Intra-defocal image. (default: {None})
            extraImg {[ndarray]} -- Extra-defocal image. (default: {None})
        """
        
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

if __name__ == "__main__":
    pass