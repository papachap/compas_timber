import math

import compas
from compas.geometry import Plane
from compas.geometry import Frame
from compas.geometry import Line
from compas.geometry import Brep
from compas.geometry import BrepTrimmingError
from compas.geometry import intersection_line_plane
from compas.geometry import distance_point_point
from compas.geometry import angle_vectors_signed
from compas.geometry import is_point_behind_plane
from compas.geometry import Vector
from compas.geometry import Rotation

from compas.tolerance import TOL

if not compas.IPY:
    from typing import TYPE_CHECKING

    if TYPE_CHECKING:
        from compas_timber.elements import Beam


class BTLxProcess(object):
    def __init__(self, ref_side_index):
        self.ref_side_index = ref_side_index


class OrientationType(object):
    START = 0
    END = 1


class JackRafterCut(BTLxProcess):
    def __init__(self, orientation, start_x=0.0, start_y=0.0, start_depth=0.0, angle=90.0, inclination=90.0, **kwargs):
        super(JackRafterCut, self).__init__(**kwargs)
        self._orientation = None
        self._start_x = None
        self._start_y = None
        self._start_depth = None
        self._angle = None
        self._inclination = None

        self.orientation = orientation
        self.start_x = start_x
        self.start_y = start_y
        self.start_depth = start_depth
        self.angle = angle
        self.inclination = inclination

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
        if start_x > 100000.0:
            raise ValueError("Start X must be less than 50000.0.")
        self._start_x = start_x

    @property
    def start_y(self):
        return self._start_y

    @start_y.setter
    def start_y(self, start_y):
        if start_y > 50000.0:
            raise ValueError("Start Y must be less than 50000.0.")
        self._start_y = start_y

    @property
    def start_depth(self):
        return self._start_depth

    @start_depth.setter
    def start_depth(self, start_depth):
        if start_depth > 50000.0:
            raise ValueError("Start Depth must be less than 50000.0.")
        self._start_depth = start_depth

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

    @classmethod
    def from_plane_and_beam(cls, plane, beam, ref_side_index=0):
        # type: (Plane, Beam, int) -> JackRafterCut
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

    def apply(self, beam, geometry):
        # type: (Beam, Brep) -> Brep
        cutting_plane = self.plane_from_params(beam)
        try:
            return geometry.trimmed(cutting_plane, TOL.absolute)
        except BrepTrimmingError:
            raise Exception(
                # TODO: use FeatureApplicationError
                "The cutting plane does not intersect with beam geometry.",
            )

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

    def plane_from_params(self, beam):
        # type: (Beam) -> Plane
        # calculates the cutting plane from the machining parameters and the beam
        # plane origin is ref
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
