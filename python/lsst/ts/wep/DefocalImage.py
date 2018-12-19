class DefocalImage(object):

    def __init__(self, intraImg=None, extraImg=None):
        """Initialize the DefocalImage class.

        Parameters
        ----------
        intraImg : numpy.ndarray, optional
            Intra-defocal image. (the default is None.)
        extraImg : numpy.ndarray, optional
            Extra-defocal image. (the default is None.)
        """

        # Defocal images
        self.intraImg = intraImg
        self.extraImg = extraImg

    def setImg(self, intraImg=None, extraImg=None):
        """Set the image.

        Parameters
        ----------
        intraImg : numpy.ndarray, optional
            Intra-defocal image. (the default is None.)
        extraImg : numpy.ndarray, optional
            Extra-defocal image. (the default is None.)
        """

        self._setValIfNotNone("intraImg", intraImg)
        self._setValIfNotNone("extraImg", extraImg)

    def _setValIfNotNone(self, attrName, val):
        """Set the value to the related class attribute if the value is not
        none.

        Parameters
        ----------
        attrName : str
            Attribute name.
        val : numpy.ndarray
            Assigned value.
        """

        if (val is not None):
            setattr(self, attrName, val)

    def getIntraImg(self):
        """Get the intra-defocal image.

        Returns
        -------
        numpy.ndarray
            Intra-defocal image.
        """

        return self.intraImg

    def getExtraImg(self):
        """Get the extra-defocal image.

        Returns
        -------
        numpy.ndarray
            Extra-defocal image.
        """

        return self.extraImg


if __name__ == "__main__":
    pass
