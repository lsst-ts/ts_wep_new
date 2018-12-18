class SliceFacade(object):

    def __init__(self, dataButlerObject):

        self.dataButlerObject = dataButlerObject

    def getDimensions(self, key=None):
        """

        Get the dimension of exposure image.

        Arguments:
            key {[string]} -- Name of amplifier. (default: {None})

        Returns:
            [float] -- Dimension of exposure image.
        """

        return self.__getItem(key).getDimensions()

    def getBoresightRotAngle(self, key=None):
        """

        Get rotation angle at boresight at middle of exposure. The meaning of rotation
        angle depends on rotType. For example, if rotType is SKY the angle is the position
        angle of the focal plane +Y with respect to North.

        Arguments:
            key {[string]} -- Name of amplifier. (default: {None})

        Returns:
            [float] -- Rotation angle at boresight.
        """

        return self.__getItem(key).getInfo().getVisitInfo().getBoresightRotAngle()

    def getRotType(self, key=None):
        """

        Get rotation type of boresightRotAngle

        Arguments:
            key {[string]} -- Name of amplifier. (default: {None})

        Returns:
            [string] -- Get the rotation type.
        """

        return self.__getItem(key).getInfo().getVisitInfo().getRotType()

    def getBoresightAzAlt(self, key=None):
        """

        Get the boresight of telescope in Az and Alt.

        Arguments:
            key {[string]} -- Name of amplifier. (default: {None})

        Returns:
            [float] -- Boresight of telescope in Az and Alt.
        """

        return self.__getItem(key).getInfo().getVisitInfo().getBoresightAzAlt()

    def getBoresightRaDec(self, key=None):
        """

        Get the boresight of telescope in Ra and Dec.

        Arguments:
            key {[string]} -- Name of amplifier. (default: {None})

        Returns:
            [float] -- Boresight of telescope in Ra and Dec.
        """

        return self.__getItem(key).getInfo().getVisitInfo().getBoresightRaDec()

    def getGain(self, key=None):
        """

        Get the gain value of amplifier.

        Arguments:
            key {[string]} -- Name of amplifier. (default: {None})

        Returns:
            [float] -- Gain of amplifier.
        """

        return self.__getItem(key).getMetadata().get("GAIN")

    def __getItem(self, key=None):
        """
        
        Get the specific exposure.
        
        Keyword Arguments:
            key {[string]} -- Name of amplifier. (default: {None})
        
        Returns:
            [butler] -- Exposure butler object.
        """

        if key is not None:
            return self.dataButlerObject[key]
        else:
            return self.dataButlerObject


if __name__ == "__main__":
    pass
