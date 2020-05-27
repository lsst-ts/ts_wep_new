# This file is part of ts_wep.
#
# Developed for the LSST Telescope and Site Systems.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import unittest

from lsst.ts.wep.Utility import CamType
from lsst.ts.wep.ctrlIntf.WEPCalculationFactory import WEPCalculationFactory
from lsst.ts.wep.ctrlIntf.WEPCalculationOfLsstCam import WEPCalculationOfLsstCam
from lsst.ts.wep.ctrlIntf.WEPCalculationOfLsstFamCam import WEPCalculationOfLsstFamCam
from lsst.ts.wep.ctrlIntf.WEPCalculationOfComCam import WEPCalculationOfComCam


class TestWEPCalculationFactory(unittest.TestCase):
    """Test the WEPCalculationFactory class."""

    def testGetCalculatorOfLsstCam(self):

        calculator = WEPCalculationFactory.getCalculator(CamType.LsstCam, "")

        self.assertTrue(isinstance(calculator, WEPCalculationOfLsstCam))

    @unittest.skip("Not support in cwfs module yet.")
    def testGetCalculatorOfLsstFamCam(self):

        calculator = WEPCalculationFactory.getCalculator(CamType.LsstFamCam, "")

        self.assertTrue(isinstance(calculator, WEPCalculationOfLsstFamCam))

    def testGetCalculatorOfComCam(self):

        calculator = WEPCalculationFactory.getCalculator(CamType.ComCam, "")

        self.assertTrue(isinstance(calculator, WEPCalculationOfComCam))

    def testGetCalculatorOfWrongCamType(self):

        self.assertRaises(
            ValueError, WEPCalculationFactory.getCalculator, "", "wrongInst"
        )


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
