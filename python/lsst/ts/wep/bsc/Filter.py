from numpy import nan
import unittest

class Filter(object):

    # LSST camera filter type
    FilterU = "u" 
    FilterG = "g"
    FilterR = "r"
    FilterI = "i"
    FilterZ = "z"
    FilterY = "y"

    def __init__(self):

        self.filter = None

    def getFilter(self):
        """
        
        Get the filter type.
        
        Returns:
            [string] -- Filter type.
        """

        return self.filter

    def setFilter(self, afilter):
        """
        
        Set the filter type.
        
        Arguments:
            afilter {[string]} -- Filter type.
        
        Raises:
            ValueError -- No such filter type.
        """

        if afilter in (self.FilterU, self.FilterG, self.FilterR, self.FilterI, self.FilterZ, self.FilterY):
            self.filter = afilter
        else:
            raise ValueError("No '%s' filter." % afilter)

    def getMagBoundary(self):
        """
        
        Get the boundary of magnitude under the current filter type.
        
        Returns:
            [float] -- Boundary of magnitude (lowMagnitude, highMagnitude).
        """

        # Get the boundary of magnitude based on the filter
        lowMagnitude = nan
        highMagnitude = nan
        if (self.filter == self.FilterU):
            lowMagnitude = 7.94
            highMagnitude = 14.80

        elif (self.filter == self.FilterG):
            lowMagnitude = 9.74
            highMagnitude = 16.17

        elif (self.filter == self.FilterR):
            lowMagnitude = 9.56
            highMagnitude = 15.73

        elif (self.filter == self.FilterI):
            lowMagnitude = 9.22
            highMagnitude = 15.26

        elif (self.filter == self.FilterZ):
            lowMagnitude = 8.83
            highMagnitude = 14.68
            
        elif (self.filter == self.FilterY):
            lowMagnitude = 8.02
            highMagnitude = 13.76

        return lowMagnitude, highMagnitude

class FilterTest(unittest.TestCase):
    """
    Test the function of Filter.
    """

    def setUp(self):

        self.filter = "u"

    def testFilter(self):

        afilter = Filter()
        afilter.setFilter(self.filter)

        self.assertEqual(afilter.getFilter(), self.filter)
        self.assertEqual(afilter.getMagBoundary(), (7.94, 14.8))

if __name__ == '__main__':

    # Do the unit test
    unittest.main()

