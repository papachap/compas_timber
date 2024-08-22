import math

from compas.geometry import Box
from compas.geometry import Brep
from compas.geometry import Frame
from compas.geometry import Line
from compas.geometry import Point
from compas.geometry import Plane
from compas.geometry import PlanarSurface
from compas.geometry import Rotation
from compas.geometry import Vector
from compas.geometry import angle_vectors_signed
from compas.geometry import distance_point_point
from compas.geometry import intersection_line_plane
from compas.geometry import is_point_behind_plane
from compas.geometry import offset_polyline
from compas.tolerance import TOL

from compas_timber.elements import FeatureApplicationError

from .btlx_process import BTLxProcess
from .btlx_process import BTLxProcessParams
from .btlx_process import OrientationType
from .btlx_process import TenonShapeType


class House(BTLxProcess):
    """Represents a House feature to be made on a beam. It should be combines with a tenon or a dovetail tenon feature.

    Parameters
    ----------
    orientation : int
        The orientation of the cut. Must be either OrientationType.START or OrientationType.END.
    start_x : float
        The start x-coordinate of the cut in parametric space of the reference side. Distance from the beam start to the reference point. -100000.0 < start_x < 100000.0.
    start_y : float
        The start y-coordinate of the cut in parametric space of the reference side. Distance from the reference edge to the reference point. -5000.0 < start_y < 5000.0.
    start_depth : float
        The start depth of the cut in parametric space of the reference side. Margin on the reference side. -5000.0 < start_depth < 5000.0.
    angle : float
        The angle of the cut. Angle between edge and reference edge. 0.1 < angle < 179.9.
    inclination : float
        The inclination of the cut. Inclination between face and reference side. 0.1 < inclination < 179.9.
    rotation : float
        The rotation of the cut. Angle between axis of the tenon and rederence side. 0.1 < rotation < 179.9.
    length_limited_top : bool
        Whether the top length of the cut is limited. True or False.
    length_limited_bottom : bool
        Whether the bottom length of the cut is limited. True or False.
    length : float
        The length of the cut. 0.0 < length < 5000.0.
    width : float
        The width of the cut. 0.0 < width < 1000.0.
    height : float
        The height of the tenon. 0.0 < height < 1000.0.
    shape : str
        The shape of the cut. Must be either 'automatic', 'square', 'round', 'rounded', or 'radius'.
    shape_radius : float
        The radius of the shape of the cut. 0.0 < shape_radius < 1000.0.

    """

    PROCESS_NAME = "House"  # type: ignore

    # Class-level attribute
    _dovetail_tool_params = {}

    @property
    def __data__(self):
        data = super(House, self).__data__
        data["orientation"] = self.orientation
        data["start_x"] = self.start_x
        data["start_y"] = self.start_y
        data["start_depth"] = self.start_depth
        data["angle"] = self.angle
        data["inclination"] = self.inclination
        data["rotation"] = self.rotation
        data["length_limited_top"] = self.length_limited_top
        data["length_limited_bottom"] = self.length_limited_bottom
        data["length"] = self.length
        data["width"] = self.width
        data["height"] = self.height
        data["shape"] = self.shape
        data["shape_radius"] = self.shape_radius
        data["chamfer"] = self.chamfer
        return data

    def __init__(
        self,
        orientation,
        start_x=0.0,
        start_y=50.0,
        start_depth=50.0,
        angle=95.0,
        inclination=10.0,
        rotation=90.0,
        length_limited_top=True,
        length_limited_bottom=True,
        length=80.0,
        width=40.0,
        height=28.0,
        shape=TenonShapeType.AUTOMATIC,
        shape_radius=20.0,
        chamfer=False,
        **kwargs
    ):
        super(House, self).__init__(**kwargs)
        self._orientation = None
        self._start_x = None
        self._start_y = None
        self._start_depth = None
        self._angle = None
        self._inclination = None
        self._rotation = None
        self._length_limited_top = None
        self._length_limited_bottom = None
        self._length = None
        self._width = None
        self._height = None
        self._shape = None
        self._shape_radius = None
        self._chamfer = None

        self.orientation = orientation
        self.start_x = start_x
        self.start_y = start_y
        self.start_depth = start_depth
        self.angle = angle
        self.inclination = inclination
        self.rotation = rotation
        self.length_limited_top = length_limited_top
        self.length_limited_bottom = length_limited_bottom
        self.length = length
        self.width = width
        self.height = height
        self.shape = shape
        self.shape_radius = shape_radius
        self.chamfer = chamfer

    ########################################################################
    # Properties
    ########################################################################

    @property
    def params_dict(self):
        return HouseParams(self).as_dict()

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
            raise ValueError("StartX must be between -100000.0 and 100000.")
        self._start_x = start_x

    @property
    def start_y(self):
        return self._start_y

    @start_y.setter
    def start_y(self, start_y):
        if start_y > 5000.0 or start_y < -5000.0:
            raise ValueError("StartY must be between -5000.0 and 5000.")
        self._start_y = start_y

    @property
    def start_depth(self):
        return self._start_depth

    @start_depth.setter
    def start_depth(self, start_depth):
        if start_depth > 5000.0 or start_depth < -5000.0:
            raise ValueError("StartDepth must be between -5000.0 and 5000.")
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

    @property
    def rotation(self):
        return self._rotation

    @rotation.setter
    def rotation(self, rotation):
        if rotation > 179.9 or rotation < 0.1:
            raise ValueError("Rotation must be between 0.1 and 179.9.")
        self._rotation = rotation

    @property
    def length_limited_top(self):
        return self._length_limited_top

    @length_limited_top.setter
    def length_limited_top(self, length_limited_top):
        if not isinstance(length_limited_top, bool):
            raise ValueError("LengthLimitedTop must be either True or False.")
        self._length_limited_top = length_limited_top

    @property
    def length_limited_bottom(self):
        return self._length_limited_bottom

    @length_limited_bottom.setter
    def length_limited_bottom(self, length_limited_bottom):
        if not isinstance(length_limited_bottom, bool):
            raise ValueError("LengthLimitedBottom must be either True or False.")
        self._length_limited_bottom = length_limited_bottom

    @property
    def length(self):
        return self._length

    @length.setter
    def length(self, length):
        if length > 5000.0 or length < 0.0:
            raise ValueError("Length must be between 0.0 and 5000.")
        self._length = length

    @property
    def width(self):
        return self._width

    @width.setter
    def width(self, width):
        if width > 1000.0 or width < 0.0:
            raise ValueError("Width must be between 0.0 and 1000.")
        self._width = width

    @property
    def height(self):
        return self._height

    @height.setter
    def height(self, height):
        if height > 1000.0 or height < 0.0:
            raise ValueError("Height must be between 0.0 and 1000.")
        self._height = height

    @property
    def shape(self):
        return self._shape

    @shape.setter
    def shape(self, shape):
        if shape not in [
            TenonShapeType.AUTOMATIC,
            TenonShapeType.SQUARE,
            TenonShapeType.ROUND,
            TenonShapeType.ROUNDED,
            TenonShapeType.RADIUS,
        ]:
            raise ValueError("Shape must be either 'automatic', 'square', 'round', 'rounded', or 'radius'.")
        self._shape = shape

    @property
    def shape_radius(self):
        return self._shape_radius

    @shape_radius.setter
    def shape_radius(self, shape_radius):
        if shape_radius > 1000.0 or shape_radius < 0.0:
            raise ValueError("ShapeRadius must be between 0.0 and 1000.")
        self._shape_radius = shape_radius

    @property
    def chamfer(self):
        return self._chamfer

    @chamfer.setter
    def chamfer(self, chamfer):
        if not isinstance(chamfer, bool):
            raise ValueError("Chamfer must be either True or False.")
        self._chamfer = chamfer

    ########################################################################
    # Alternative constructors
    ########################################################################

    @classmethod
    def from_plane_and_beam(
        cls,
        plane,
        beam,
        start_depth=50.0,
        rotation=90.0,
        length=80.0,
        width=40.0,
        height=28.0,
        shape=TenonShapeType.AUTOMATIC,
        shape_radius=20.0,
        chamfer=False,
        ref_side_index=0,
    ):
        """Create a House instance from a cutting surface and the beam it should cut. This could be the ref_side of the cross beam of a Joint and the main beam.

        Parameters
        ----------
        plane : :class:`~compas.geometry.Plane` or :class:`~compas.geometry.Frame`
            The cutting plane.
        beam : :class:`~compas_timber.elements.Beam`
            The beam that is cut by this instance.
        ref_side_index : int, optional
            The reference side index of the beam to be cut. Default is 0 (i.e. RS1).

        Returns
        -------
        :class:`~compas_timber.fabrication.House`

        """
        # type: (Plane|Frame, Beam, float, float, bool, int) -> House

        if isinstance(plane, Frame):
            plane = Plane.from_frame(plane)
        # define ref_side & ref_edge
        ref_side = beam.ref_sides[ref_side_index]
        ref_edge = Line.from_point_and_vector(ref_side.point, ref_side.xaxis)

        # calculate orientation
        orientation = cls._calculate_orientation(ref_side, plane)

        # calculate start_x
        point_start_x = intersection_line_plane(ref_edge, plane)
        if point_start_x is None:
            raise ValueError("Plane does not intersect with beam.")
        start_x = distance_point_point(ref_side.point, point_start_x)

        # calculate start_y
        start_y = beam.width / 2  # TODO: is there a case where the tenon should not be in the middle of the beam?

        # calculate angle
        angle = cls._calculate_angle(ref_side, plane, orientation)

        # calculate inclination
        inclination = cls._calculate_inclination(ref_side, plane, orientation)

        # calculate rotation
        rotation = 90.0  # TODO: for now i'm overriding unit a logic for applying rotation is implemented

        # bound start_depth, length and width
        start_depth = cls._bound_start_depth(orientation, start_depth, inclination, height)
        length = cls._bound_length(beam.height, start_depth, inclination, length, height)
        width = cls._bound_width(beam.width, width, shape_radius)

        # determine if the top and bottom length of the cut is limited
        # TODO: this should instead come first and override the start_depth and length
        length_limited_top = start_depth > 0.0  #!
        length_limited_bottom = length < ((beam.height) / math.sin(math.radians(inclination)) - start_depth)  #!
        # override because otherwise tenon would go out of the blank
        if inclination > 90.0:
            length_limited_bottom = True
            # length_limited_bottom = length < max_length - frustum_difference
        elif inclination < 90.0:
            length_limited_top = True

        return cls(
            orientation,
            start_x,
            start_y,
            start_depth,
            angle,
            inclination,
            rotation,
            length_limited_top,
            length_limited_bottom,
            length,
            width,
            height,
            shape,
            shape_radius,
            chamfer,
            ref_side_index=ref_side_index,
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

    @staticmethod
    # bound the start_depth value to the minimum possible start_depth if the incliantion is larger than 90
    def _bound_start_depth(orientation, start_depth, inclination, height):
        if orientation == OrientationType.START:
            min_start_depth = height / (math.tan(math.radians(inclination)))
        else:
            min_start_depth = height / (math.tan(math.radians(180 - inclination)))
        return max(start_depth, min_start_depth)

    @staticmethod
    def _bound_length(beam_height, start_depth, inclination, length, height):
        # bound the inserted lenhgth value to the maximum possible length for the beam based on the inclination
        max_length = (beam_height) / math.sin(math.radians(inclination)) - start_depth
        print("max_length", max_length)
        if inclination < 90.0:
            max_length = max_length - height / math.tan(math.radians(inclination))
        return min(max_length, length)

    @staticmethod
    def _bound_width(beam_width, width, shape_radius):
        # bound the inserted width value to the minumum and maximum possible width for the beam
        min_width = 2 * shape_radius
        return min(max(width, min_width), beam_width)

    @staticmethod
    def _create_trimming_planes_for_box(bottom_points, top_points, length_limited_top, length_limited_bottom):
        # Create the trimming planes to trim a box into a frustum trapezoid.
        if len(bottom_points) != len(top_points):
            raise ValueError("The number of bottom points must match the number of top points.")

        planes = []
        num_points = len(bottom_points)
        for i in range(num_points):
            # Get the bottom and top points for the current edge
            bottom1 = bottom_points[i]
            top1 = top_points[i]
            bottom2 = bottom_points[(i + 1) % num_points]
            # Define vectors for the plane
            x_axis = top1 - bottom1
            y_axis = bottom2 - bottom1
            # Compute the normal vector to the plane using the cross product
            normal = x_axis.cross(y_axis)
            # Create the plane
            plane = Plane(bottom1, -normal)
            planes.append(plane)

        return planes

    ########################################################################
    # Methods
    ########################################################################

    def apply(self, geometry, beam):
        """Apply the feature to the beam geometry.

        Parameters
        ----------
        geometry : :class:`compas.geometry.Brep`
            The geometry to be processed.

        beam : :class:`compas_timber.elements.Beam`
            The beam that is milled by this instance.

        Raises
        ------
        :class:`~compas_timber.elements.FeatureApplicationError`
            If the cutting planes do not create a volume that itersects with beam geometry or any step fails.

        Returns
        -------
        :class:`~compas.geometry.Brep`
            The resulting geometry after processing

        """
        # type: (Brep, Beam) -> Brep

        # get cutting planes from params and beam
        try:
            cutting_plane = self.plane_from_params_and_beam(beam)
            cutting_plane.normal = -cutting_plane.normal
        except ValueError as e:
            raise FeatureApplicationError(
                None, geometry, "Failed to generate cutting plane from parameters and beam: {}".format(str(e))
            )

        # get house volume from params and beam
        try:
            house_volume = self.house_volume_from_params_and_beam(beam)
        except ValueError as e:
            raise FeatureApplicationError(
                None, geometry, "Failed to generate house volume from parameters and beam: {}".format(str(e))
            )

        # fillet the edges of the house volume based on the shape
        if (
            self.shape not in [TenonShapeType.SQUARE, TenonShapeType.AUTOMATIC] and not self.length_limited_bottom
        ):  # TODO: Remove AUTOMATIC once Brep Fillet is implemented
            edge_ideces = [4, 7] if self.length_limited_top else [5, 8]
            try:
                dovetail_volume.fillet(
                    self.shape_radius, [dovetail_volume.edges[edge_ideces[0]], dovetail_volume.edges[edge_ideces[1]]]
                )  # TODO: NotImplementedError
            except Exception as e:
                raise FeatureApplicationError(
                    dovetail_volume,
                    geometry,
                    "Failed to fillet the edges of the dovetail volume based on the shape: {}".format(str(e)),
                )

        # trim geometry with cutting planes
        try:
            geometry.trim(cutting_plane)
        except Exception as e:
            raise FeatureApplicationError(
                cutting_plane, geometry, "Failed to trim geometry with cutting plane: {}".format(str(e))
            )

        # add tenon volume to geometry
        try:
            geometry += dovetail_volume
        except Exception as e:
            raise FeatureApplicationError(
                dovetail_volume, geometry, "Failed to add tenon volume to geometry: {}".format(str(e))
            )

        return geometry

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

    def house_volume_from_params_and_beam(self, beam):
        """Calculates the house volume from the machining parameters in this instance and the given beam.

        Parameters
        ----------
        beam : :class:`compas_timber.elements.Beam`
            The beam that is cut by this instance.

        Returns
        -------
        :class:`compas.geometry.Brep`
            The tenon volume.

        """
        # type: (Beam) -> Brep

        assert self.inclination is not None
        assert self.height is not None
        assert self.shape is not None
        assert self.shape_radius is not None
        assert self.length_limited_top is not None
        assert self.length_limited_bottom is not None

        cutting_frame = Frame.from_plane(self.plane_from_params_and_beam(beam))
        if self.orientation == OrientationType.START:
            cutting_frame.xaxis = -cutting_frame.xaxis
        cutting_surface = PlanarSurface(
            beam.height / math.sin(math.radians(self.inclination)), beam.width, cutting_frame
        )

        dx_top = self.length * math.tan(math.radians(self.cone_angle))
        dx_bottom = (beam.width - self.width) / 2

        bottom_house_points = [
            cutting_surface.point_at(self.start_y - dx_top, -self.start_depth),
            cutting_surface.point_at(self.start_y + dx_top, -self.start_depth),
            cutting_surface.point_at(dx_bottom + self.width, -self.start_depth - self.length),
            cutting_surface.point_at(dx_bottom, -self.start_depth - self.length),
            cutting_surface.point_at(self.start_y - dx_top, -self.start_depth),
        ]

        # Calculate the offset length
        offset_length = self.height * math.tan(math.radians(self.flank_angle))
        # offset the polyline to create the top face of the tenon
        top_house_points = offset_polyline(bottom_dovetail_points, offset_length, cutting_frame.normal)
        # make the top face flat
        top_dovetail_points[0][2] = bottom_dovetail_points[0][2]
        top_dovetail_points[1][2] = bottom_dovetail_points[1][2]
        # remove the last point to avoid duplication
        top_dovetail_points = [Point(pt[0], pt[1] + self.height, pt[2]) for pt in top_dovetail_points[:-1]]
        bottom_dovetail_points.pop(-1)

        # create the house volume by trimming a box  # TODO: PluginNotInstalledError for Brep.from_loft
        # get the box as a brep
        frame_for_box = cutting_frame.copy()
        frame_for_box.point = cutting_surface.point_at(
            beam.width / 2, -(beam.height / 2) / math.sin(math.radians(self.inclination))
        )
        dovetail_volume = Brep.from_box(
            Box(beam.width, beam.height / math.sin(math.radians(self.inclination)) * 1.5, self.height, frame_for_box)
        )

        translation_vector = cutting_frame.normal * (self.height / 2)
        house_volume.translate(translation_vector)

        # get trimming planes for creating the dovetail volume
        trimming_planes = self._create_trimming_planes_for_box(
            bottom_dovetail_points, top_dovetail_points, self.length_limited_top, self.length_limited_bottom
        )

        # trim the box to create the dovetail volume
        for plane in trimming_planes:
            try:
                if self.orientation == OrientationType.START:
                    plane.normal = -plane.normal
                dovetail_volume.trim(plane)
            except Exception as e:
                raise FeatureApplicationError(
                    plane, dovetail_volume, "Failed to trim tenon volume with cutting plane: {}".format(str(e))
                )

        return dovetail_volume


class HouseParams(BTLxProcessParams):
    """A class to store the parameters of a house feature.

    Parameters
    ----------
    instance : :class:`~compas_timber._fabrication.House`
        The instance of the house feature.
    """

    def __init__(self, instance):
        # type: (House) -> None
        super(HouseParams, self).__init__(instance)

    def as_dict(self):
        """Returns the parameters of the House feature as a dictionary.

        Returns
        -------
        dict
            The parameters of the House feature as a dictionary.
        """
        # type: () -> OrderedDict
        result = super(HouseParams, self).as_dict()
        result["Orientation"] = self._instance.orientation
        result["StartX"] = "{:.{prec}f}".format(self._instance.start_x, prec=TOL.precision)
        result["StartY"] = "{:.{prec}f}".format(self._instance.start_y, prec=TOL.precision)
        result["StartDepth"] = "{:.{prec}f}".format(self._instance.start_depth, prec=TOL.precision)
        result["Angle"] = "{:.{prec}f}".format(self._instance.angle, prec=TOL.precision)
        result["Inclination"] = "{:.{prec}f}".format(self._instance.inclination, prec=TOL.precision)
        result["Rotation"] = "{:.{prec}f}".format(self._instance.rotation, prec=TOL.precision)
        result["LengthLimitedTop"] = "yes" if self._instance.length_limited_top else "no"
        result["LengthLimitedBottom"] = "yes" if self._instance.length_limited_bottom else "no"
        result["Length"] = "{:.{prec}f}".format(self._instance.length, prec=TOL.precision)
        result["Width"] = "{:.{prec}f}".format(self._instance.width, prec=TOL.precision)
        result["Height"] = "{:.{prec}f}".format(self._instance.height, prec=TOL.precision)
        result["ConeAngle"] = "{:.{prec}f}".format(self._instance.cone_angle, prec=TOL.precision)
        result["UseFlankAngle"] = "yes" if self._instance.use_flank_angle else "no"
        result["FlankAngle"] = "{:.{prec}f}".format(self._instance.flank_angle, prec=TOL.precision)
        result["Shape"] = self._instance.shape
        result["ShapeRadius"] = "{:.{prec}f}".format(self._instance.shape_radius, prec=TOL.precision)
        return result
