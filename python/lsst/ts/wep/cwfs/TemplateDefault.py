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


class TemplateDefault(object):
    """Default Template class."""

    def makeTemplate(self, sensorName, defocalState, imageSize, **kwargs):
        """Make the template image.

        Parameters
        ----------
        sensorName : str
            The camera detector for which we want to make a template. Should
            be in "Rxx_Sxx" format.
        defocalState : str
            "extra" or "intra" describing the defocal state of the sensor.
        imageSize : int
            Size of template in pixels. The template will be a square.
        **kwargs : dict[str, any]
            Dictionary of input argument: new value for that input argument.

        Returns
        -------
        NotImplementedError
            Child class should implement this.
        """

        raise NotImplementedError("Child class should implement this.")
