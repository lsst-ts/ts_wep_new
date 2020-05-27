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

from lsst.ts.wep.ctrlIntf.WEPCalculation import WEPCalculation


class WEPCalculationOfPiston(WEPCalculation):
    """The child class of WEPCalculation that gets the defocal images by the
    camera piston."""

    def __init__(self, astWcsSol, camType, isrDir):
        """Construct an WEP calculation of piston object.

        Parameters
        ----------
        astWcsSol : AstWcsSol
            AST world coordinate system (WCS) solution.
        camType : enum 'CamType'
            Camera type.
        isrDir : str
            Instrument signature remocal (ISR) directory. This directory will
            have the input and output that the data butler needs.
        """
        super(WEPCalculationOfPiston, self).__init__(astWcsSol, camType, isrDir)

    def setDefocalDisInMm(self, defocalDisInMm):
        """Set the defocal distance in mm.

        Parameters
        ----------
        float
            Defocal distance in mm.
        """

        wfEsti = self.wepCntlr.getWfEsti()

        inst = wfEsti.getInst()
        inst.setAnnDefocalDisInMm(defocalDisInMm)

    def getDefocalDisInMm(self):
        """Set the defocal distance in mm.

        CWFS: Curvature wavefront sensor.

        Returns
        -------
        float
            Defocal distance in mm used in the cwfs algorithm.
        """

        wfEsti = self.wepCntlr.getWfEsti()

        inst = wfEsti.getInst()
        defocalDisInM = inst.getDefocalDisOffset()

        return defocalDisInM * 1e3


if __name__ == "__main__":
    pass
