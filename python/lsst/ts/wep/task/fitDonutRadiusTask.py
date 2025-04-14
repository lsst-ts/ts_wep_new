import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import numpy as np
from astropy.table import QTable
from lsst.pipe.base import connectionTypes as ct
from lsst.ts.wep.task.donutStamps import DonutStamps
from lsst.utils.timer import timeMethod
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks

__all__ = [
    "FitDonutRadiusTaskConnections",
    "FitDonutRadiusTaskConfig",
    "FitDonutRadiusTask",
]


class FitDonutRadiusTaskConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("visit", "detector", "instrument"),
):

    donutStampsIntra = ct.Input(
        doc="Intra-focal Donut Stamps",
        dimensions=("visit", "detector", "instrument"),
        storageClass="StampsBase",
        name="donutStampsIntra",
    )
    donutStampsExtra = ct.Input(
        doc="Extra-focal Donut Stamps",
        dimensions=("visit", "detector", "instrument"),
        storageClass="StampsBase",
        name="donutStampsExtra",
    )
    donutRadiiTable = ct.Output(
        doc="Donut radius table",
        dimensions=("visit", "detector", "instrument"),
        storageClass="AstropyQTable",
        name="donutRadiiTable",
    )


class FitDonutRadiusTaskConfig(
    pipeBase.PipelineTaskConfig,
    pipelineConnections=FitDonutRadiusTaskConnections,
):

    widthMultiplier = pexConfig.Field[float](
        doc="Multiplier used to convert the width of peaks fitted \
         to donut cross-section to donut edge.",
        default=0.8,
    )
    filterSigma = pexConfig.Field[float](
        doc="Standard deviation of the Gaussian kernel \
        used to smooth out the donut cross-section prior to \
        using peak finder (in pixels).",
        default=3,
    )
    minPeakWidth = pexConfig.Field[float](
        doc="Required minimum width of peaks (in pixels) in \
        donut cross-section.",
        default=5,
    )
    minPeakHeight = pexConfig.Field[float](
        doc="Required minimum height of peaks in normalized \
        donut cross-section (i.e. must be between 0 and 1).",
        default=0.3,
    )
    leftDefault = pexConfig.Field[float](
        doc="Default position of the left edge of the donut \
        expressed as a fraction of image length (i.e. between \
        0 and 1).",
        default=0.1,
    )
    rightDefault = pexConfig.Field[float](
        doc="Default position of the right edge of the donut \
        expressed as a fraction of image length (i.e. between \
        0 and 1).",
        default=0.8,
    )


class FitDonutRadiusTask(pipeBase.PipelineTask):
    ConfigClass = FitDonutRadiusTaskConfig
    _DefaultName = "FitDonutRadius"
    """
    Run donut radius fitting.
    """

    def empty(self) -> pipeBase.Struct:
        """Return empty table if no donut stamps are available.

        Parameters
        ----------
        donutRadiiTable : astropy.table.QTable
            Donut radius table.

        Returns
        -------
        lsst.pipe.base.Struct
            Empty donut radii table.
        """

        donutRadiiTable = QTable(
            {
                "VISIT": [],
                "DFC_TYPE": [],
                "DET_NAME": [],
                "DFC_DIST": [],
                "RADIUS": [],
                "X_LEFT_EDGE": [],
                "X_RIGHT_EDGE": [],
                "FAIL_FLAG": [],
            }
        )

        return pipeBase.Struct(donutRadiiTable=donutRadiiTable)

    @timeMethod
    def run(
        self,
        donutStampsExtra: DonutStamps,
        donutStampsIntra: DonutStamps,
    ) -> pipeBase.Struct:
        # If no donuts are in the DonutStamps
        # Then return an empty table
        if len(donutStampsExtra) == 0 or len(donutStampsIntra) == 0:
            return self.empty()

        widthMultiplier = self.config.widthMultiplier
        filterSigma = self.config.filterSigma
        minPeakWidth = self.config.minPeakWidth
        minPeakHeight = self.config.minPeakHeight
        leftDefault = self.config.leftDefault
        rightDefault = self.config.rightDefault

        dfcTypeArray = []
        visitArray = []
        radiiArray = []
        leftEdgeArray = []
        rightEdgeArray = []
        centX0Array = []
        centY0Array = []
        detNameArray = []
        defocalDistArray = []
        failedFlagArray = []

        # dfcType = stamps.metadata["DFC_TYPE"]  # single string
        # stampVisit = stamps.metadata["VISIT"]  # single int

        # looping over each  defocal pair
        for stampSet in [donutStampsIntra, donutStampsExtra]:
            dfcType = stampSet.metadata["DFC_TYPE"]  # single string
            stampVisit = stampSet.metadata["VISIT"]  # single int

            dfcTypeArray.append(np.full(len(stampSet), dfcType, dtype="<U5"))
            visitArray.append(np.ones(len(stampSet)) * stampVisit)
            centX0Array.append(stampSet.metadata.getArray("CENT_X0"))  # list of values
            centY0Array.append(stampSet.metadata.getArray("CENT_Y0"))  # list of values

            for stamp in stampSet:
                image = stamp.stamp_im.image.array
                xLeft, xRight, donutRadius, flag = self.fit_radius(
                    image=image,
                    multiplier=widthMultiplier,
                    filter_sigma=filterSigma,
                    min_peak_width=minPeakWidth,
                    min_height=minPeakHeight,
                    left_default_perc=leftDefault,
                    right_default_perc=rightDefault,
                )
                radiiArray.append(donutRadius)
                leftEdgeArray.append(xLeft)
                rightEdgeArray.append(xRight)
                detNameArray.append(stamp.detector_name)
                defocalDistArray.append(stamp.defocal_distance)
                failedFlagArray.append(flag)

        dfcTypeArrayFlat = np.concatenate([np.asarray(item) for item in dfcTypeArray])
        visitArrayFlat = np.concatenate(
            [np.asarray(item) for item in visitArray]
        ).astype(int)

        donutRadiiTable = QTable(
            {
                "VISIT": visitArrayFlat,
                "DFC_TYPE": dfcTypeArrayFlat,
                "DET_NAME": detNameArray,
                "DFC_DIST": defocalDistArray,
                "RADIUS": radiiArray,
                "X_LEFT_EDGE": leftEdgeArray,
                "X_RIGHT_EDGE": rightEdgeArray,
                "FAIL_FLAG": failedFlagArray,
            }
        )

        return pipeBase.Struct(donutRadiiTable=donutRadiiTable)

    def fit_radius(
        self,
        image,
        multiplier,
        filter_sigma,
        min_peak_width,
        min_height,
        left_default_perc,
        right_default_perc,
    ):
        """
        Find peaks and widths of smoothed donut cross-section.
        The donut cross-section is normalized, and
        smoothed out with Gaussian kernel
        prior to peak finding. The location and width
        of each peak is used to determine the left
        and right edge of the donut.

        Parameters
        ----------
        image : numpy.ndarray
            The donut stamp image.
        multiplier: float
            Multiplier used to convert the width of peaks fitted
            to donut cross-section to donut edge.
        filter_sigma: float
            Standard deviation of the Gaussian kernel
            used to smooth out the donut cross-section prior to
            using peak finder (in pixels).",
        min_peak_width: float
            Required minimum width of peaks (in pixels) in
            donut cross-section.
        min_height: float
            Required minimum height of peaks in normalized
            donut cross-section (i.e. must be between 0 and 1).
        left_default_perc: float
            The position of the left edge of the donut in case
            fitting failed: a fraction of image length
            (between 0 and 1).
        right_default_perc: float
            The position of the right edge of the donut in case
            fitting failed: a fraction of image length
            (between 0 and 1).

        Returns
        -------
        left_edge: float
            The position of the left donut edge in pixels.
        right_edge: float
            The position of the right donut edge in pixels.
        radius: float
            The donut radius in pixels.
        fail_flag: int
            Flag set to 1 if any fit failure encountered
            (in which case default values are stored, and
            a warning is logged).

        """
        halfWidth = int(len(image) / 2)
        y_cross = image[halfWidth, :]
        y_cross_norm = np.array(y_cross) / max(y_cross)
        fail_flag = 0
        # convolve with Gaussian filter to smooth out the cross-section
        filtered_x_profile = gaussian_filter(y_cross_norm, sigma=filter_sigma)

        # set the defaults used in case of fit failure,
        # so that the returned radius will be reasonable
        left_default_edge = left_default_perc * len(image)
        right_default_edge = right_default_perc * len(image)

        # detect sources
        peak_locations, peak_information = find_peaks(
            filtered_x_profile,
            height=min_height,
            width=min_peak_width,
        )
        if len(peak_locations) > 0:

            # choose left and right peaks
            index_of_right = np.argmax(peak_locations)
            index_of_left = np.argmin(peak_locations)

            left_width = peak_information["widths"][index_of_left]
            right_width = peak_information["widths"][index_of_right]
            left_peak = peak_locations[index_of_left]
            right_peak = peak_locations[index_of_right]

            left_edge = left_peak - left_width * multiplier
            right_edge = right_peak + right_width * multiplier
        else:
            left_edge = left_default_edge
            right_edge = right_default_edge
            fail_flag = 1
            self.log.warning(
                f"Setting left edge to {left_edge} and right edge to {right_edge}"
            )

        # Catch successful fit with bad values
        if left_edge < 0:
            self.log.warning(f"warning:  left_edge: {left_edge}")
            left_edge = left_default_edge
            fail_flag = 1
        if right_edge > len(image):
            self.log.warning(f"warning:  right_edge: {right_edge}")
            right_edge = right_default_edge
            fail_flag = 1

        # donut radius is half of the distance
        # between two edges
        radius = (right_edge - left_edge) / 2.0
        return left_edge, right_edge, radius, fail_flag
