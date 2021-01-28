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

import os
import numpy as np
from lsst.ts.wep.Utility import getConfigDir, readPhoSimSettingData, CamType
from lsst.ts.wep.cwfs.TemplateDefault import TemplateDefault
from lsst.ts.wep.cwfs.Instrument import Instrument
from lsst.ts.wep.cwfs.CompensableImage import CompensableImage


class TemplateModel(TemplateDefault):
    """Class to make the donut templates from the Instrument model."""
    def makeTemplate(self, sensorName, defocalState, imageSize,
                     camType=CamType.LsstCam, pixelScale=0.2):
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
        camType : enum 'CamType'
            Camera type. (Default is CamType.LsstCam)
        pixelScale : float
            The pixels to arcseconds conversion factor. (The default is 0.2)

        Returns
        -------
        numpy.ndarray
            The donut template as a binary image
        """

        configDir = getConfigDir()
        focalPlaneLayout = readPhoSimSettingData(configDir, 'focalplanelayout.txt', "fieldCenter")

        pixelSizeInUm = float(focalPlaneLayout[sensorName][2])
        sizeXinPixel = int(focalPlaneLayout[sensorName][3])

        sensorXMicron, sensorYMicron = np.array(focalPlaneLayout[sensorName][:2], dtype=float)
        # Correction for wavefront sensors
        # (from _shiftCenterWfs in SourceProcessor.py)
        if sensorName in ("R44_S00_C0", "R00_S22_C1"):
            # Shift center to +x direction
            sensorXMicron = sensorXMicron + sizeXinPixel / 2 * pixelSizeInUm
        elif sensorName in ("R44_S00_C1", "R00_S22_C0"):
            # Shift center to -x direction
            sensorXMicron = sensorXMicron - sizeXinPixel / 2 * pixelSizeInUm
        elif sensorName in ("R04_S20_C1", "R40_S02_C0"):
            # Shift center to -y direction
            sensorYMicron = sensorYMicron - sizeXinPixel / 2 * pixelSizeInUm
        elif sensorName in ("R04_S20_C0", "R40_S02_C1"):
            # Shift center to +y direction
            sensorYMicron = sensorYMicron + sizeXinPixel / 2 * pixelSizeInUm

        # Load Instrument parameters
        instDir = os.path.join(configDir, "cwfs", "instData")
        dimOfDonutImgOnSensor = imageSize
        inst = Instrument(instDir)
        inst.config(camType, dimOfDonutImgOnSensor)

        # Create image for mask
        img = CompensableImage()
        img.defocalType = defocalState

        # Convert pixel locations to degrees
        sensorXPixel = float(sensorXMicron)/pixelSizeInUm
        sensorYPixel = float(sensorYMicron)/pixelSizeInUm

        sensorXDeg = sensorXPixel*pixelScale / 3600
        sensorYDeg = sensorYPixel*pixelScale / 3600

        # define position of donut at center of current sensor in degrees
        boundaryT = 0
        maskScalingFactorLocal = 1
        img.fieldX, img.fieldY = sensorXDeg, sensorYDeg
        img.makeMask(inst, "offAxis", boundaryT, maskScalingFactorLocal)

        templateArray = img.cMask

        return templateArray
