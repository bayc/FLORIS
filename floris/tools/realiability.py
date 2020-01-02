import numpy as np

class Vec3():
    def __init__(self, x1, x2=None, x3=None, string_format=None):
        """
        Object containing vector information for coordinates.

        Args:
            x1: [float, float, float] or float -- The first argument
            can be a list of the three vector components or simply the
            first component of the vector.
            x2: float (optional) -- The second component of the vector.
            x3: float (optional) -- The third component of the vector.
            string_format: str (optional) -- The string format to use in the
                overloaded __str__ function.
        """
        if isinstance(x1, list):
            self.x1, self.x2, self.x3 = [x for x in x1]
        else:
            self.x1 = x1
            self.x2 = x2
            self.x3 = x3

        if not (type(self.x1) == type(self.x2) and
                type(self.x1) == type(self.x3) and
                type(self.x2) == type(self.x3)):
                target_type = type(self.x1)
                self.x2 = target_type(self.x2)
                self.x3 = target_type(self.x3)

        if string_format is not None:
            self.string_format = string_format
        else:
            if type(self.x1) in [int]:
                self.string_format = "{:8d}"
            elif type(self.x1) in [float, np.float64]:
                self.string_format = "{:8.3f}"

    def rotate_on_x3(self, theta, center_of_rotation=None):
        """
        Rotate about the x3 coordinate axis by a given angle
        and center of rotation.
        The angle theta should be given in degrees.

        Sets the rotated components on this object and returns
        """
        if center_of_rotation is None:
            center_of_rotation = Vec3(0.0, 0.0, 0.0)
        x1offset = self.x1 - center_of_rotation.x1
        x2offset = self.x2 - center_of_rotation.x2
        self.x1prime = x1offset * cosd(theta) - x2offset * sind(
            theta) + center_of_rotation.x1
        self.x2prime = x2offset * cosd(theta) + x1offset * sind(
            theta) + center_of_rotation.x2
        self.x3prime = self.x3

    def __str__(self):
        template_string = "{} {} {}".format(self.string_format,
                                            self.string_format,
                                            self.string_format)
        return template_string.format(self.x1, self.x2, self.x3)

    def __add__(self, arg):

        if type(arg) is Vec3:
            return Vec3(self.x1 + arg.x1, self.x2 + arg.x2, self.x3 + arg.x3)
        else:
            return Vec3(self.x1 + arg, self.x2 + arg, self.x3 + arg)

    def __sub__(self, arg):
        if type(arg) is Vec3:
            return Vec3(self.x1 - arg.x1, self.x2 - arg.x2, self.x3 - arg.x3)
        else:
            return Vec3(self.x1 - arg, self.x2 - arg, self.x3 - arg)

    def __mul__(self, arg):

        if type(arg) is Vec3:
            return Vec3(self.x1 * arg.x1, self.x2 * arg.x2, self.x3 * arg.x3)
        else:
            return Vec3(self.x1 * arg, self.x2 * arg, self.x3 * arg)

    def __truediv__(self, arg):
        if type(arg) is Vec3:
            return Vec3(self.x1 / arg.x1, self.x2 / arg.x2, self.x3 / arg.x3)
        else:
            return Vec3(self.x1 / arg, self.x2 / arg, self.x3 / arg)

    def __eq__(self, arg):
        return self.x1 == arg.x1 \
            and self.x2 == arg.x2 \
            and self.x3 == arg.x3

    def __hash__(self):
        return hash((self.x1, self.x2, self.x3))


def get_surface(wind_speed, TI, T1, T2, component):
    """
    This function gets the respective surrogate surface for the input 
    parameters.
    
    Args:
        wind_speed ([type]): [description]
        TI ([type]): [description]
        T1 ([type]): [description]
        T2 ([type]): [description]
        component ([type]): [description]
    """

    # pickle.load(the_right_file)

# this would live in turbine.py
def calculate_turbine_reliablity(wind_speed, TI, T1, T2, component):
    surf = get_surface(wind_speed, TI, T1, T2, component)

    reliability = calc_reliability(x, y, surf, component)

    self.component_reliability = reliability

# this would run after calculate_wake, using the flow field 

from itertools import combinations

def calculate_farm_reliability(wind_speed, TI, component, fi, surf_dir):

    with open(surf_dir + 'Case_A_8_L10_T1', 'rb') as dill_file:
        L10_T1_surface = pickle.load(dill_file)
    
    with open(surf_dir + 'Case_A_8_L10_T2', 'rb') as dill_file:
        L10_T2_surface = pickle.load(dill_file)

    # Turbine_L10s = np.zeros(self.nturbs)
    Turbine_L10s = np.zeros(4)
    # Turbine_PowerVals = np.zeros(self.nturbs)

    # turbs = []
    # theta = math.radians(wind_direction) #rotate grid 90 deg
    # for ii in list(zip(layout_x, layout_y)):
    #     turbs.append((self.rotate_array(ii, theta)))
    
    # layout_x = [ii[0] for ii in turbs]
    # layout_y = [ii[1] for ii in turbs]

    # rotated_map = self.turbine_map.rotated(
    #         self.wind_direction, center_of_rotation)

    center_of_rotation = Vec3(0, 0, 0)

    rotated_map = fi.floris.farm.flow_field.turbine_map.rotated(
            0, center_of_rotation)

    # sort the turbine map
    sorted_map = rotated_map.sorted_in_x_as_list()

    layout_x = [coord.x1 for coord, turbine in sorted_map]
    layout_y = [coord.x2 for coord, turbine in sorted_map]
    
    turb_label = range(0, len(layout_x))    

    # 2) Assign L10 values to all turbines
    
    print(L10_T1_surface[0])
    L10_matrix = np.full( (len(turb_label), len(turb_label)),
                           float(L10_T1_surface[0]) )
    # Power_matrix = np.full( (len(turb_label), len(turb_label)),
    #                          float(Power_T1_surface[0]) )

    turbs_to_analyse = []

    x_threshold = 1260
    y_threshold = 252
    
    #Delete pairs that fall outside space constraints
    for ii in list(combinations(turb_label, 2)):
        turb1 = ii[0] 
        turb2 = ii[1]
        if abs(layout_x[turb2] - layout_x[turb1]) <= x_threshold and \
           abs(layout_y[turb2] - layout_y[turb1]) <= y_threshold:
            turbs_to_analyse.append(ii)
            
    # Assign L10 and power values; turbs are ordered in proximity to the
    # incoming wind direction

    for ii in turbs_to_analyse:
        turb1 = ii[0]
        turb2 = ii[1]
        L10_matrix[turb1][turb2] = L10_T2_surface(layout_y[turb2] \
                                   - layout_y[turb1], layout_x[turb2] \
                                   - layout_x[turb1])
        # Power_matrix[turb1][turb2] = Power_T2_surface(layout_x[turb2],
        #                                               layout_y[turb2]) #kW/hr

    #Find the minimum of each row (this is the minimum )
    for column in L10_matrix.T:
        Turbine_L10s = [min(column) for column in L10_matrix.T]
        # Turbine_PowerVals = [Power_matrix.T[ii][np.argmin(column)] \
        #                      for ii, column in enumerate(L10_matrix.T)]

    return Turbine_L10s


import floris.tools as wfct
import dill as pickle

fi = wfct.floris_utilities.FlorisInterface('../../examples/example_input.json')

surf_dir = '../../examples/'

wind_speed = []
TI = []
component = []

tmp = calculate_farm_reliability(wind_speed, TI, component, fi, surf_dir)

print(tmp)