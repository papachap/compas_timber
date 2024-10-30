import math

from compas.geometry import BrepTrimmingError
from compas.geometry import Frame
from compas.geometry import Line
from compas.geometry import Plane
from compas.geometry import Rotation
from compas.geometry import Vector
from compas.geometry import angle_vectors_signed
from compas.geometry import distance_point_point
from compas.geometry import intersection_line_plane
from compas.geometry import is_point_behind_plane
from compas.tolerance import TOL

from compas_timber.elements import FeatureApplicationError

from .btlx_process import BTLxProcess
from .btlx_process import BTLxProcessParams
from .btlx_process import OrientationType


class Lap(BTLxProcess):
    """Represents a Lap feature to be made on a beam.

    Parameters
    ----------
    orientation : int
        The orientation of the cut. Must be either OrientationType.START or OrientationType.END.
    start_x : float
        The start x-coordinate of the cut in parametric space of the reference side. -100000.0 < start_x < 100000.0.
    start_y : float
        The start y-coordinate of the cut in parametric space of the reference side. -50000.0 < start_y < 50000.0.
    angle : float
        The horizontal angle of the cut. 0.1 < angle < 179.9.
    inclination : float
        The vertical angle of the cut. 0.1 < inclination < 179.9.
    slope : float
        The slope of the cut. -89.9 < slope < 89.9.
    length : float
        The length of the cut. 0.0 < length < 100000.0.
    width : float
        The width of the cut. 0.0 < width < 50000.0.
    depth : float
        The depth of the cut. -50000.0 < depth < 50000.0.
    lead_angle_parallel : bool
        The lead angle is parallel to the beam axis.
    lead_angle : float
        The lead angle of the cut. 0.1 < lead_angle < 179.9.
    lead_inclination_parallel : bool
        The lead inclination is parallel to the beam axis.
    lead_inclination : float
        The lead inclination of the cut. 0.1 < lead_inclination < 179.9.
    machining_limits : dict, optional
        The machining limits for the cut. Default is None

    """

    PROCESS_NAME = "Lap"  # type: ignore

    @property
    def __data__(self):
        data = super(Lap, self).__data__
        data["orientation"] = self.orientation
        data["start_x"] = self.start_x
        data["start_y"] = self.start_y
        data["angle"] = self.angle
        data["inclination"] = self.inclination
        data["slope"] = self.slope
        data["length"] = self.length
        data["width"] = self.width
        data["depth"] = self.depth
        data["lead_angle_parallel"] = self.lead_angle_parallel
        data["lead_angle"] = self.lead_angle
        data["lead_inclination_parallel"] = self.lead_inclination_parallel
        data["lead_inclination"] = self.lead_inclination
        data["machining_limits"] = self.machining_limits
        return data

    # fmt: off
    def __init__(
        self,
        orientation,
        start_x=0.0,
        start_y=0.0,
        angle=90.0,
        inclination=90.0,
        slope=0.0,
        length=200.0,
        width=50.0,
        depth=40.0,
        lead_angle_parallel=True,
        lead_angle=90.0,
        lead_inclination_parallel=True,
        lead_inclination=90.0,
        machining_limits=None,
        **kwargs
    ):
        super(Lap, self).__init__(**kwargs)
        self._orientation = None
        self._start_x = None
        self._start_y = None
        self._angle = None
        self._inclination = None
        self._slope = None
        self._length = None
        self._width = None
        self._depth = None
        self._lead_angle_parallel = None
        self._lead_angle = None
        self._lead_inclination_parallel = None
        self._lead_inclination = None
        self._machining_limits = None

        self.orientation = orientation
        self.start_x = start_x
        self.start_y = start_y
        self.angle = angle
        self.inclination = inclination
        self.slope = slope
        self.length = length
        self.width = width
        self.depth = depth
        self.lead_angle_parallel = lead_angle_parallel
        self.lead_angle = lead_angle
        self.lead_inclination_parallel = lead_inclination_parallel
        self.lead_inclination = lead_inclination
        self.machining_limits = machining_limits

    ########################################################################
    # Properties
    ########################################################################

    @property
    def params_dict(self):
        return LapParams(self).as_dict()

    @property
    def orientation(self):
        return self._orientation

    @orientation.setter
    def orientation(self, orientation):
        if orientation not in [OrientationType.START, OrientationType.END]:
            raise ValueError("Orientation must be either OrientationType.START or OrientationType.END.")
        self._orientation = orientation

    @property
    def start_x(self):
        return self._start_x

    @start_x.setter
    def start_x(self, start_x):
        if start_x > 100000.0 or start_x < -100000.0:
            raise ValueError("Start X must be between -100000.0 and 100000.")
        self._start_x = start_x

    @property
    def start_y(self):
        return self._start_y

    @start_y.setter
    def start_y(self, start_y):
        if -50000.0 > start_y > 50000.0:
            raise ValueError("Start Y must be between -50000.0 and 50000.")
        self._start_y = start_y

    @property
    def angle(self):
        return self._angle

    @angle.setter
    def angle(self, angle):
        if angle > 179.9 or angle < 0.1:
            raise ValueError("Angle must be between 0.1 and 179.9.")
        self._angle = angle

    @property
    def inclination(self):
        return self._inclination

    @inclination.setter
    def inclination(self, inclination):
        if inclination > 179.9 or inclination < 0.1:
            raise ValueError("Inclination must be between 0.1 and 179.9.")
        self._inclination = inclination

    @property
    def slope(self):
        return self._slope

    @slope.setter
    def slope(self, slope):
        if slope > 89.9 or slope < -89.9:
            raise ValueError("Slope must be between -89.9 and 89.9.")
        self._slope = slope

    @property
    def length(self):
        return self._length

    @length.setter
    def length(self, length):
        if length > 100000.0 or length < 0.0:
            raise ValueError("Length must be between 0.0 and 100000.")
        self._length = length

    @property
    def width(self):
        return self._width

    @width.setter
    def width(self, width):
        if width > 50000.0 or width < 0.0:
            raise ValueError("Width must be between 0.0 and 50000.")
        self._width = width

    @property
    def depth(self):
        return self._depth

    @depth.setter
    def depth(self, depth):
        if depth > 50000.0 or depth < -50000.0:
            raise ValueError("Depth must be between -50000 and 50000.")
        self._depth = depth

    @property
    def lead_angle_parallel(self):
        return self._lead_angle_parallel

    @lead_angle_parallel.setter
    def lead_angle_parallel(self, lead_angle_parallel):
        if not isinstance(lead_angle_parallel, bool):
            raise ValueError("Lead angle parallel must be a boolean.")
        self._lead_angle_parallel = lead_angle_parallel

    @property
    def lead_angle(self):
        return self._lead_angle

    @lead_angle.setter
    def lead_angle(self, lead_angle):
        if lead_angle > 179.9 or lead_angle < 0.1:
            raise ValueError("Lead angle must be between 0.1 and 179.9.")
        self._lead_angle = lead_angle

    @property
    def lead_inclination_parallel(self):
        return self._lead_inclination_parallel

    @lead_inclination_parallel.setter
    def lead_inclination_parallel(self, lead_inclination_parallel):
        if not isinstance(lead_inclination_parallel, bool):
            raise ValueError("Lead inclination parallel must be a boolean.")
        self._lead_inclination_parallel = lead_inclination_parallel

    @property
    def lead_inclination(self):
        return self._lead_inclination

    @lead_inclination.setter
    def lead_inclination(self, lead_inclination):
        if lead_inclination > 179.9 or lead_inclination < 0.1:
            raise ValueError("Lead inclination must be between 0.1 and 179.9.")
        self._lead_inclination = lead_inclination

    @property
    def machining_limits(self):
        return self._machining_limits

    @machining_limits.setter
    def machining_limits(self, machining_limits):
        if machining_limits is not None and not isinstance(machining_limits, dict):
            raise ValueError("Machining limits must be a dictionary.")
        self._machining_limits = machining_limits


    ########################################################################
    # Alternative constructors
    ########################################################################

    @classmethod
    def from_planes_and_beam(cls, plane, beam, depth, ref_side_index=0):
        """Create a Lap instance from a set of cutting planes and the beam it should cut.

        Parameters
        ----------
        planes : list of :class:`~compas.geometry.Plane` or :class:`~compas.geometry.Frame`
            The cutting planes to be used for the Lap feature.
        beam : :class:`~compas_timber.elements.Beam`
            The beam that is cut by this instance.
        depth : float
            The depth of the lap.
        ref_side_index : int, optional
            The reference side index of the beam to be cut. Default is 0 (i.e. RS1).

        Returns
        -------
        :class:`~compas_timber.fabrication.Lap`

        """
        # type: (list[Plane | Frame], Beam, float, int) -> Lap

        if isinstance(plane, Frame):
            plane = Plane.from_frame(plane)
        start_y = 0.0
        start_depth = 0.0
        ref_side = beam.ref_sides[ref_side_index]  # TODO: is this arbitrary?
        ref_edge = Line.from_point_and_vector(ref_side.point, ref_side.xaxis)
        orientation = cls._calculate_orientation(ref_side, plane)

        point_start_x = intersection_line_plane(ref_edge, plane)
        if point_start_x is None:
            raise ValueError("Plane does not intersect with beam.")

        start_x = distance_point_point(ref_edge.point, point_start_x)
        angle = cls._calculate_angle(ref_side, plane, orientation)
        inclination = cls._calculate_inclination(ref_side, plane, orientation)
        return cls(orientation, start_x, start_y, start_depth, angle, inclination, ref_side_index=ref_side_index)

    @staticmethod
    def _calculate_orientation(ref_side, cutting_plane):
        # orientation is START if cutting plane normal points towards the start of the beam and END otherwise
        # essentially if the start is being cut or the end
        if is_point_behind_plane(ref_side.point, cutting_plane):
            return OrientationType.END
        else:
            return OrientationType.START

    @staticmethod
    def _calculate_angle(ref_side, plane, orientation):
        # vector rotation direction of the plane's normal in the vertical direction
        angle_vector = Vector.cross(ref_side.zaxis, plane.normal)
        angle = angle_vectors_signed(ref_side.xaxis, angle_vector, ref_side.zaxis, deg=True)
        if orientation == OrientationType.START:
            return 180 - abs(angle)  # get the other side of the angle
        else:
            return abs(angle)

    @staticmethod
    def _calculate_inclination(ref_side, plane, orientation):
        # vector rotation direction of the plane's normal in the horizontal direction
        inclination_vector = Vector.cross(ref_side.yaxis, plane.normal)
        inclination = angle_vectors_signed(ref_side.xaxis, inclination_vector, ref_side.yaxis, deg=True)
        if orientation == OrientationType.START:
            return 180 - abs(inclination)  # get the other side of the angle
        else:
            return abs(inclination)

    ########################################################################
    # Methods
    ########################################################################

    def apply(self, geometry, beam):
        """Apply the feature to the beam geometry.

        Parameters
        ----------
        geometry : :class:`~compas.geometry.Brep`
            The beam geometry to be cut.
        beam : :class:`compas_timber.elements.Beam`
            The beam that is cut by this instance.

        Raises
        ------
        :class:`~compas_timber.elements.FeatureApplicationError`
            If the cutting plane does not intersect with beam geometry.

        Returns
        -------
        :class:`~compas.geometry.Brep`
            The resulting geometry after processing

        """
        # type: (Brep, Beam) -> Brep
        cutting_plane = self.plane_from_params_and_beam(beam)
        try:
            return geometry.trimmed(cutting_plane)
        except BrepTrimmingError:
            raise FeatureApplicationError(
                cutting_plane,
                geometry,
                "The cutting plane does not intersect with beam geometry.",
            )

    def plane_from_params_and_beam(self, beam):
        """Calculates the cutting plane from the machining parameters in this instance and the given beam

        Parameters
        ----------
        beam : :class:`compas_timber.elements.Beam`
            The beam that is cut by this instance.

        Returns
        -------
        :class:`compas.geometry.Plane`
            The cutting plane.

        """
        # type: (Beam) -> Plane
        assert self.angle is not None
        assert self.inclination is not None

        # start with a plane aligned with the ref side but shifted to the start_x of the cut
        ref_side = beam.side_as_surface(self.ref_side_index)
        p_origin = ref_side.point_at(self.start_x, 0.0)
        cutting_plane = Frame(p_origin, ref_side.frame.xaxis, ref_side.frame.yaxis)

        # normal pointing towards xaxis so just need the delta
        horizontal_angle = math.radians(self.angle - 90)
        rot_a = Rotation.from_axis_and_angle(cutting_plane.zaxis, horizontal_angle, point=p_origin)

        # normal pointing towards xaxis so just need the delta
        vertical_angle = math.radians(self.inclination - 90)
        rot_b = Rotation.from_axis_and_angle(cutting_plane.yaxis, vertical_angle, point=p_origin)

        cutting_plane.transform(rot_a * rot_b)
        # for simplicity, we always start with normal pointing towards xaxis.
        # if start is cut, we need to flip the normal
        if self.orientation == OrientationType.END:
            plane_normal = cutting_plane.xaxis
        else:
            plane_normal = -cutting_plane.xaxis
        return Plane(cutting_plane.point, plane_normal)


class LapParams(BTLxProcessParams):
    """A class to store the parameters of a Lap feature.

    Parameters
    ----------
    instance : :class:`~compas_timber._fabrication.Lap`
        The instance of the Lap feature.
    """

    def __init__(self, instance):
        # type: (Lap) -> None
        super(LapParams, self).__init__(instance)

    def as_dict(self):
        """Returns the parameters of the Lap feature as a dictionary.

        Returns
        -------
        dict
            The parameters of the Lap feature as a dictionary.
        """
        # type: () -> OrderedDict
        result = super(LapParams, self).as_dict()
        result["Orientation"] = self._instance.orientation
        result["StartX"] = "{:.{prec}f}".format(self._instance.start_x, prec=TOL.precision)
        result["StartY"] = "{:.{prec}f}".format(self._instance.start_y, prec=TOL.precision)
        result["StartDepth"] = "{:.{prec}f}".format(self._instance.start_depth, prec=TOL.precision)
        result["Angle"] = "{:.{prec}f}".format(self._instance.angle, prec=TOL.precision)
        result["Inclination"] = "{:.{prec}f}".format(self._instance.inclination, prec=TOL.precision)
        return result
