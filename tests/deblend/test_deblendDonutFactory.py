import unittest

from lsst.ts.wep.deblend.DeblendDonutFactory import DeblendDonutFactory
from lsst.ts.wep.deblend.DeblendAdapt import DeblendAdapt
from lsst.ts.wep.Utility import DeblendDonutType


class TestDeblendDonutFactory(unittest.TestCase):
    """Test the DeblendDonutFactory class."""

    def testCreateDeblendAdapt(self):

        deblendDonut = DeblendDonutFactory.createDeblendDonut(
            DeblendDonutType.Adapt)
        self.assertTrue(isinstance(deblendDonut, DeblendAdapt))

    def testCreateDeblendDonutWrongType(self):

        self.assertRaises(ValueError, DeblendDonutFactory.createDeblendDonut,
                          "wrongType")


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
