import math
import numpy as np
import openmdao.api as om
#-------------------------------------------------------------------------------
def read_data_input(file):
    """
    DESCRIPTION: Reads the specified inputs from the data input file and stores
                 the information in a Python dictionary for its easy access
                 later on when executing code routines and subroutines.
    INPUTS:
        * file: File from which the input data is read ---> "Data_Input.txt"
    OUTPUTS:
        * {data_dictionary}: Python dictionary containing the read and stored
                           input data. Data is stored by pairs of key-item. The
                           key to access every item respect the same name that
                           the input from the "Data_Input.txt". For further
                           information about Python dictionaries the user is
                           encouraged to consult Section 5.5: Data structures:
                           Dictionaries from the Python documentation, available
                           in the link:
                           https://docs.python.org/3/tutorial/datastructures.html#dictionaries

    """
    infile = open(file, "r")
    data_input = infile.read().split()
    infile.close()

    data_dictionary = {
                       "total_forces_file_name"     : data_input[64],
                       "shear_moment_file_name"     : data_input[66],
                       "strip_forces_file_name"     : data_input[68],
                       "airfoil_points_file_name"   : data_input[71],
                       "altitude"                   : float(data_input[74]),
                       "velocity"                   : float(data_input[76]),
                       "reference_surface"          : float(data_input[78]),
                       "span"                       : float(data_input[80]),
                       "stringer_yield_stress"      : float(data_input[83]),
                       "skin_yield_stress"          : float(data_input[85]),
                       "stringer_elastic_modulus"   : float(data_input[87]),
                       "skin_elastic_modulus"       : float(data_input[89]),
                       "stringer_material_density"  : float(data_input[91]),
                       "skin_material_density"      : float(data_input[93]),
                       "stringer_type"              : (data_input[96]).upper(),
                       "web_minimum_length"         : float(data_input[98]),
                       "web_maximum_length"         : float(data_input[100]),
                       "flange_length_proportion"   : float(data_input[102]),
                       "web_minimum_thickness"      : float(data_input[104]),
                       "web_maximum_thickness"      : float(data_input[106]),
                       "flange_thickness_proportion": float(data_input[108]),
                       "stringer_reserve_factor"    : float(data_input[110]),
                       "skin_minimum_thickness"     : float(data_input[112]),
                       "skin_maximum_thickness"     : float(data_input[114]),
                       "skin_reserve_factor"        : float(data_input[116])
                       }

    return data_dictionary
#-------------------------------------------------------------------------------
def extract_chords(data_dictionary):
    """
    DESCRIPTION: To reduce the amount of inputs required to the user, chords at
                 root and tip position are extracted from the AVL file
                 "strip_forces_file_name" which contains some geometrical data.
    INPUTS:
        * {data_dictionary}: The dictionary containing all the user input data.
                           All filenames and variables can be obtained from here.
    OUTPUTS: There is no output produces, there is no 'return' in the function.
             The 'data_dictionary' is updated with two key-item pairs containing
             the value of chord at root and at tip: data_dictionary["chord_root"]
             and data_dictionary["chord_tip"].
    """
    file = data_dictionary["strip_forces_file_name"]
    infile = open(file, "r")
    data_lines = infile.readlines()
    infile.close()

    for line in range(len(data_lines)):
        if data_lines[line] == "\n":
            continue
        elif data_lines[line].split()[0] == "j":
            init = line + 1
        elif data_lines[line].split()[0] == "Surface" and \
             data_lines[line].split()[1] == "#" and \
             data_lines[line].split()[2] == "2":
             end = line - 2
             break

    data_dictionary["chord_root"] = float(data_lines[init].split()[2])
    data_dictionary["chord_tip"] = float(data_lines[end].split()[2])
#-------------------------------------------------------------------------------
def standard_atmosphere(altitude):
    """
    DESCRIPTION: From the altitude and using the ISA model (International
                 Standard Atmosphere) calculates the air density, necessary
                 for obtaining the dynamic pressure and extract the real loads
                 and aerodynamic forces from their non-dimensional form.
                 The function includes a Troposphere model valid until 11.000
                 meters, a lower Stratosphere model valid from 11.000 to
                 25.000 meters and an upper Stratosphere model valid for altitudes
                 higher than 25.000 meters. For the further information about the
                 implementation of the models, the user is encouraged to consult
                 the ISA model by NASA, available in the link:
                 https://www.grc.nasa.gov/WWW/K-12/airplane/atmosmet.html
    INPUTS:
        * altitude: Flight altitude in meters.
    OUTPUTS:
        * rho: Air density at the corresponding flight altitude.
    """
    if altitude <= 11000:
        T = 15.04 - 0.00649*altitude #Temperature in ºC
        p = 101.29 * ((T+273.1)/288.08)**5.256 #Pressure in KPa
        rho = p / (0.2869*(T+273.1)) #Density in kg/m^3

    elif altitude > 11000 and altitude <= 25000:
        T = -56.46 #Temperature in ºC
        p = 22.65 * exp(1.73-0.000157*altitude) #Pressure in KPa
        rho = p / (0.2869*(T+273.1)) #Density in kg/m^3

    else:
        T = -131.21 + 0.00299*altitude #Temperature in ºC
        p = 2.488 * ((T+273.1)/216.6)**-11.388 #Pressure in KPa
        rho = p / (0.2869*(T+273.1)) #Density in kg/m^3

    return rho
#-------------------------------------------------------------------------------
def extract_shear_moment(data_dictionary):
    """
    DESCRIPTION: Extracts the list of shear force and bending moment variation
                 along the wingspan from the AVL file "shear_moment_file_name".
                 AVL data is non-dimensionalised, data is dimensionalised
                 using the user input data about flight conditions.
                 Data is stored in lists for later access during calculations in
                 code routines and subroutines.
    INPUTS:
        * {data_dictionary}: The dictionary containing all the user input data.
                           All filenames and variables can be obtained from here.
    OUTPUTS:
        * [Y1_list]: A Python list containing the discrete span positions.
        * [Vz_list]: A Python list containing the vertical shear force at the
                   span positions.
        * [Mx_list]: A Python list containing the bending moment at the span
                   positions.
    """
    file = data_dictionary["shear_moment_file_name"]
    infile = open(file, "r")
    data_lines = infile.readlines()
    infile.close()

    for line in range(len(data_lines)):
        if data_lines[line] == "\n":
            continue
        elif data_lines[line].split()[0] == "2Y/Bref":
            init = line + 1
        elif data_lines[line] == " Surface:   2\n":
            end = line - 2
            break
    h = data_dictionary["altitude"]
    V = data_dictionary["velocity"]
    Sref = data_dictionary["reference_surface"]
    Bref = data_dictionary["span"]
    rho = standard_atmosphere(h)
    q = (1/2)*rho*V**2 #Dynamic pressure

    Y1_list = []
    Vz_list = []
    Mx_list = []
    for line in range(init,end+1):
        Y1_list.append(float(data_lines[line].split()[0]) * (Bref/2))
        Vz_list.append(float(data_lines[line].split()[1]) * (q*Sref))
        Mx_list.append(float(data_lines[line].split()[2]) * (q*Bref*Sref))

    return Y1_list, Vz_list, Mx_list
#-------------------------------------------------------------------------------
def extract_drag(data_dictionary):
    """
    DESCRIPTION: Extracts the drag coefficient from the "total_forces_file_name".
                 Data is dimensionalised using the user input data about
                 the flight conditions and final value of drag is obtained.
    INPUTS:
        * {data_dictionary}: The dictionary containing all the user input data.
                           All filenames and variables can be obtained from here.
    OUTPUTS:
        * Vx: Total drag value.
    """
    file = data_dictionary["total_forces_file_name"]
    infile = open(file, "r")
    data_lines = infile.readlines()
    infile.close()

    for line in range(len(data_lines)):
        if data_lines[line] == "\n":
            continue
        elif data_lines[line].split()[0] == "CDtot":
            Cd = float(data_lines[line].split()[2])
            break

    h = data_dictionary["altitude"]
    V = data_dictionary["velocity"]
    Sref = data_dictionary["reference_surface"]
    Bref = data_dictionary["span"]
    rho = standard_atmosphere(h)
    q = (1/2)*rho*V**2 #Dynamic pressure
    Vx = q*Sref*Cd

    return Vx
#-------------------------------------------------------------------------------
def extract_torsion(data_dictionary):
    """
    DESCRIPTION: Extracts the list of torsion moment variation along the
                 wingspan from the AVL file "strip_forces_file_name". AVL data
                 is non-dimensionalised, data is dimensionalised using the
                 user input data about flight conditions. Data is stored
                 in lists for later access during calculations in code routines
                 and subroutines.
    INPUTS:
        * {data_dictionary}: The dictionary containing all the user input data.
                           All filenames and variables can be obtained from here.
    OUTPUTS:
        * [Y2_list]: A Python list containing the discrete span positions.
        * [My_list]: A Python list containing the torsion moment at the span
                   positions.
    """
    file = data_dictionary["strip_forces_file_name"]
    infile = open(file, "r")
    data_lines = infile.readlines()
    infile.close()

    for line in range(len(data_lines)):
        if data_lines[line] == "\n":
            continue
        elif data_lines[line].split()[0] == "j":
            init = line + 1
        elif data_lines[line].split()[0] == "Surface" and \
             data_lines[line].split()[1] == "#" and \
             data_lines[line].split()[2] == "2":
             end = line - 2
             break

    h = data_dictionary["altitude"]
    V = data_dictionary["velocity"]
    Sref = data_dictionary["reference_surface"]
    Bref = data_dictionary["span"]
    rho = standard_atmosphere(h)
    q = (1/2)*rho*V**2 #Dynamic pressure
    cref = Sref / Bref #Reference chord

    My_list = []
    Y2_list = []
    for line in range(init, end+1):
        Y2_list.append(float(data_lines[line].split()[1]))
        My_list.append(float(data_lines[line].split()[10]) * (q*Sref*cref))

    return Y2_list, My_list
#-------------------------------------------------------------------------------
def extract_airfoil(data_dictionary, scale):
    """
    DESCRIPTION: Read the cartesian coordinates of the points that conform the
                 the section under study.

    INPUTS:
        * {data_dictionary}: The dictionary containing all the user input data.
                           All filenames and variables can be obtained from here.
        * scale: The scaling factor applied to all the section, with this
                 parameter, characteristic dimensions can be adjusted and
                 units changed.
                 For the case of a wing profile, through scaling, the size of the
                 profile can be changed. If profile points are obtained such that
                 the reference chord length is 1m, the new chord length will be
                 obtained by the expression: chord = scale * 1, therefore,
                 the scale parameter will indicate the chord of the
                 airfoil.

    OUTPUTS:
        * [x]: Python list containing the x coordinates of the section.
        * [z]: Python list containing the y coordinates of the section.
    """
    file = data_dictionary["airfoil_points_file_name"]
    infile = open(file, "r")
    points = infile.read().split()
    infile.close()

    x = []
    z = []
    i = 0
    while i <= (len(points) - 1):
        x.append(scale * float(points[i]))
        z.append(scale * float(points[i + 1]))

        i += 2

    return x, z
#-------------------------------------------------------------------------------
def gravity_centre(x, z):
    """
    DESCRIPTION: Returns the position of the centre of gravity of the profile.

    INPUTS:
        * [x]: Python list containing the x coordinates of the section.
        * [z]: Python list containing the y coordinates of the section.

    OUTPUTS:
        * xg: x coordinate of the centre of gravity.
        * zg: y coordinate of the centre of gravity.
    """
    x_num = 0
    z_num = 0
    for i in range(len(x) - 1):
        x_num += x[i]
        z_num += z[i]

    xg = x_num / (len(x)-1)
    zg = z_num / (len(x)-1)

    return xg, zg
#-------------------------------------------------------------------------------
def change_to_CG(x, z, xg, zg):
    """
    DESCRIPTION: Changes the reference frame to an absolute reference frame
                 which origin of coordinates is taken in the center of gravity.
                 This allows simplifications when executing algorithms where
                 the distance between a stringer and the gravity center is required.

    INPUTS:
        * [x]: Python list containing the x coordinates of the section.
        * [z]: Python list containing the y coordinates of the section.
        * xg: x position of the gravity centre in the initial XY reference frame.
        * zg: y position of the gravity centre in the initial XY reference frame.

    OUTPUTS:
        * [x_cg]: Python list containing the x coordinates of the section with
                  respect a reference frame centred in (xg, yg).
        * [z_cg]: Python list containing the x coordinates of the section with
                  respect a reference frame centred in (xg, yg).
    """
    x_cg = []
    z_cg = []
    for i in range(len(x)):
        x_cg.append(x[i] - xg)
        z_cg.append(z[i] - zg)

    return x_cg, z_cg
#-------------------------------------------------------------------------------
def interpolation_shear_moment(y, Y1_list, Vz_list, Mx_list):
    """
    DESCRIPTION: AVL files provide a list of discrete span positions where shear
                 and bending moment are calculated. Given such lists of data,
                 this function calculates the load at any other arbitraty point
                 of the span by interpolation of the already existing data.
    INPUTS:
        * y: Span position where loads are desired to be obtained.
        * [Y1_list]: Discrete span postions provided by AVL for shear and bending
                   moment.
        * [Vz_list]: List of shear force at the the discrete span positions.
        * [Mx_list]: List of bending moment at the discrete span positions.
    OUTPUTS:
        * Vz: Shear value at the input span position.
        * Mx: Bending moment value at the input wingspan position.

    """
    exact = -1
    for n in range(len(Y1_list)):
        if y == Y1_list[n]:
            #Y1 = Y1_list[n]
            Vz = Vz_list[n]
            Mx = Mx_list[n]
            exact = n
            break

    if exact == -1:
        for n in range(1, len(Y1_list)):
            if y > Y1_list[n-1] and y < Y1_list[n]:
                Delta = Y1_list[n] - Y1_list[n-1]
                l = y - Y1_list[n-1]
                u = Y1_list[n] - y
                l_coef = 1 - l/Delta
                u_coef = 1 - u/Delta

                #Y1 = Y1_list[n-1]*l_coef + Y1_list[n]*u_coef #Y = y checking purpose
                Vz = Vz_list[n-1]*l_coef + Vz_list[n]*u_coef
                Mx = Mx_list[n-1]*l_coef + Mx_list[n]*u_coef
                break

    return Vz, Mx
#-------------------------------------------------------------------------------
def interpolation_torsion(y, Y2_list, My_list):
    """
    DESCRIPTION: AVL files provide a list of discrete span positions where
                 torsional moment is calculated. Given such lists of data,
                 this function calculates the load at any other arbitraty point
                 of the span by interpolation of the already existing data.
    INPUTS:
        * y: Span position where loads are desired to be obtained.
        * [Y2_list]: Discrete span postions provided by AVL for shear and bending
                   moment.
        * [My_list]: List of torsional moment at the discrete span positions.
    """
    exact = -1
    for n in range(len(Y2_list)):
        if y == Y2_list[n]:
            #Y2 = Y2_list[n]
            My = My_list[n]
            exact = n
            break

    if exact == -1:
        for n in range(1, len(Y2_list)):
            if y > Y2_list[n-1] and y < Y2_list[n]:
                Delta = Y2_list[n] - Y2_list[n-1]
                l = y - Y2_list[n-1]
                u = Y2_list[n] - y
                l_coef = 1 - l/Delta
                u_coef = 1 - u/Delta

                #Y1 = Y1_list[n-1]*l_coef + Y1_list[n]*u_coef #Y = y checking purpose
                My = My_list[n-1]*l_coef + My_list[n]*u_coef
                break

    return My
#-------------------------------------------------------------------------------
def stringer_area(type, web_length, flange_lenght_proportion, \
                  web_thickness, flange_thickness_proportion):
    """
    DESCRIPTION: From the stringer geometrical properties and stringer kind,
                 the function calculates the area of the stringer in order to
                 carry out the boom simplification.

    INPUTS:
        * type: Type of profile: C, Z, T or I.
        * web_length: Web length.
        * flange_lenght_proportion: Ratio between flange and web lenght.
                                    L_flange / L_web
        * web_thickness: Web thickness
        * flange_thickness_proportion: Ratio between flange and web thickness.
                                       t_flange / t_web

    OUTPUTS:
        * A_stringer: Stringer area.
    """
    flange_lenght = flange_lenght_proportion * web_length
    flange_thickness = flange_thickness_proportion * web_thickness

    if type == "C" or type == "Z" or type == "I":
        A_stringer = 2*flange_lenght*flange_thickness + web_length*web_thickness

    elif type == "T":
        A_stringer = flange_lenght*flange_thickness + web_length*web_thickness

    return A_stringer
#-------------------------------------------------------------------------------
def chord_variation_law(data_dictionary, y):
    """
    DESCRIPTION: Given a span position calculates the chord length assuming a
                 linear chord variation.
    INPUTS:
        * {data_dictionary}: The dictionary containing all the user input data.
                           All filenames and variables can be obtained from here.
        * y: Span position where loads are desired to be obtained.
    OUTPUTS:
        * c: chord length at input span position.
    """
    c_root = data_dictionary["chord_root"]
    c_tip = data_dictionary["chord_tip"]
    b = data_dictionary["span"]

    c = c_root - (c_root - c_tip)/(b/2) * y

    return c
#-------------------------------------------------------------------------------
def stringer_preliminary_sizing(data_dictionary, y, Y1_list, Vz_list, Mx_list):
    """
    DESCRIPTION: Sizing the stringer dimensions and stringer number to withstand
                 the loading at position "y" of the wingspan.
    INPUTS:
        * {data_dictionary}: The dictionary containing all the user input data.
                       All filenames and variables can be obtained from here.
        * y: Span position where loads are desired to be obtained.
        * [Y1_list]: A Python list containing the discrete span positions.
        * [Vz_list]: A Python list containing the vertical shear force at the
                   span positions.
        * [Mx_list]: A Python list containing the bending moment at the span
                   positions.
    OUTPUTS:
        * Mx: The moment at the winspan position "y". It is provided as an
              output in order for its further use and not to be interpolated again.
        * z_max: Maximum wing section height where bending effect is higher.
        * N_stringers: The minimum number of stringers needed to withsand the
                       loads caused by the bending moment.
        * web_l: Optimized web length.
        * web_t: Optimized web thickness.
        * flange_l: Optimized flange length.
        * flange_t: Optimized flange thickness.
    """
    type = data_dictionary["stringer_type"]
    web_l_min = data_dictionary["web_minimum_length"]
    web_l_max = data_dictionary["web_maximum_length"]
    flange_l_prop = data_dictionary["flange_length_proportion"]
    web_t_min = data_dictionary["web_minimum_thickness"]
    web_t_max = data_dictionary["web_maximum_thickness"]
    flange_t_prop = data_dictionary["flange_thickness_proportion"]

    Asmin = stringer_area(type, web_l_min, flange_l_prop, web_t_min, flange_t_prop)
    Asmax = stringer_area(type, web_l_max, flange_l_prop, web_t_min, flange_t_prop)

    c = chord_variation_law(data_dictionary, y)
    x, z = extract_airfoil(data_dictionary, c)
    xg, zg = gravity_centre(x, z)
    x_cg, z_cg = change_to_CG(x, z, xg, zg)
    z_max = max(z_cg)*(10**3) # To mm
    Vz, Mx = interpolation_shear_moment(y, Y1_list, Vz_list, Mx_list)
    Mx = Mx*(10**3)
    stringer_yield = data_dictionary["stringer_yield_stress"]
    stringer_RF = data_dictionary["stringer_reserve_factor"]

    class Stringer(om.ExplicitComponent):
        def setup(self):
            #Declare input design variables
            self.add_input('N', val = 0.0)
            self.add_input('As', val = 0.0)
            #Declare output variables
            self.add_output('A_tot', val = 0.0)
            # Finite difference all partials
            self.declare_partials('*', '*', method = 'fd')

        def compute(self, inputs, outputs):
            N = inputs['N']
            As = inputs['As']
            outputs['A_tot'] = N*As

    prob = om.Problem()
    indeps = prob.model.add_subsystem('indeps', om.IndepVarComp())
    indeps.add_output('N', 1)
    indeps.add_output('As', Asmin)

    prob.model.add_subsystem('Stringer_Opt', Stringer())
    prob.model.add_subsystem('Constraint', om.ExecComp('g = s_yield - RF*Mx/(N*As*z)', \
                                                        s_yield = stringer_yield, \
                                                        RF = stringer_RF, \
                                                        Mx = Mx, \
                                                        z = z_max))
    prob.model.connect('indeps.N', ['Stringer_Opt.N', 'Constraint.N'])
    prob.model.connect('indeps.As', ['Stringer_Opt.As', 'Constraint.As'])
    prob.driver = om.ScipyOptimizeDriver()
    prob.driver.options['optimizer'] = 'SLSQP'
    prob.model.add_design_var('indeps.N', lower = 1)
    prob.model.add_design_var('indeps.As', lower = Asmin, upper = Asmax)
    prob.model.add_objective('Stringer_Opt.A_tot')
    prob.model.add_constraint('Constraint.g', lower = 0)
    prob.setup()
    prob.run_driver()

    A_stringer = float(prob['indeps.As'])
    N_stringers = math.ceil(float(prob['indeps.N']))

    if A_stringer <= ((Asmax+Asmin)/2):
        factor = A_stringer/Asmin
        if type == "C" or type == "Z" or type == "I":
            web_l = factor*web_l_min
            web_t = A_stringer/(web_l*(1+2*flange_l_prop*flange_t_prop))
            flange_l = web_l*flange_l_prop
            flange_t = web_t*flange_t_prop
        elif type == "T":
            web_l = factor*web_l_min
            web_t = A_stringer/(web_l*(1+flange_l_prop*flange_t_prop))
            flange_l = web_l*flange_l_prop
            flange_t = web_t*flange_t_prop

    elif A_stringer > ((Asmax+Asmin)/2):
        factor = A_stringer/Asmax
        if type == "C" or type == "Z" or type == "I":
            web_l = factor*web_l_max
            web_t = A_stringer/(web_l*(1+2*flange_l_prop*flange_t_prop))
            flange_l = web_l*flange_l_prop
            flange_t = web_t*flange_t_prop
        elif type == "T":
            web_l = factor*web_l_max
            web_t = A_stringer/(web_l*(1+flange_l_prop*flange_t_prop))
            flange_l = web_l*flange_l_prop
            flange_t = web_t*flange_t_prop

    return Mx, z_max, N_stringers, web_l, web_t, flange_l, flange_t
#-------------------------------------------------------------------------------
def stringer_crippling(data_dictionary, Mx, z_max, N_stringers, web_length, \
                       web_thickness):
    """
    DESCRIPTION: Checks for stringer crippling given the loading conditions of
                 the section and given the preliminary optimized stringer
                 dimensions. If crippling detected changes the section dimensions
                 between the allowable values to prevent it. If crippling is not
                 possible to avoid, provides with a warning message.
    INPUTS:
        * {data_dictionary}: The dictionary containing all the user input data.
                       All filenames and variables can be obtained from here.
        * Mx: The moment at the winspan position "y". It is provided as an
              output in order for its further use and not to be interpolated again.
        * z_max: Maximum wing section height where bending effect is higher.
        * N_stringers: The minimum number of stringers needed to withsand the
                       loads caused by the bending moment.
        * web_length: Optimized web length.
        * web_thickness: Optimized web thickness.
    OUTPUTS:
        * web_length: Web length after crippling check
        * web_thickness:Web thickness after crippling check.
        * flange_length: Flange length crippling check.
        * flange_thickness: Flange thickness crippling check.
    """
    E = data_dictionary["stringer_elastic_modulus"]*(10**3)
    s_yield = data_dictionary["stringer_yield_stress"]
    flange_l_prop = data_dictionary["flange_length_proportion"]
    flange_t_prop = data_dictionary["flange_thickness_proportion"]
    B = 4.05
    m = 0.82
    type = data_dictionary["stringer_type"]
    Area_stringer = stringer_area(type, web_length, flange_l_prop, web_thickness, \
                                  flange_t_prop)
    RF = data_dictionary["stringer_reserve_factor"]

    if type == "Z" or type == "I" or type == "C":
        flange_length = web_length * flange_l_prop
        flange_thickness = web_thickness * flange_t_prop
        mean_t = web_thickness * (web_length*web_thickness)/Area_stringer + \
                + flange_thickness * (2*flange_length*flange_thickness)/Area_stringer
    elif type == "T":
        flange_length = web_length * flange_l_prop
        flange_thickness = web_thickness * flange_t_prop
        mean_t = web_thickness * (web_length*web_thickness)/Area_stringer + \
                + flange_thickness * (flange_length*flange_thickness)/Area_stringer

    load = Mx / (N_stringers*Area_stringer*z_max*1000)
    crippling_cond = ((load*RF) / (B*(E**m * s_yield**(3-m))**(1/3)))**(1/m)

    if crippling_cond <= (mean_t**2 / Area_stringer):
        print('Stringer dimensions are not prone to suffer stringer crippling.')
        print('Stringer sizing is finished.')
    elif crippling_cond > (mean_t**2 / Area_stringer):
        print('Stringer crippling detected, stringer dimensions will be ' +
              'recalculated.')

        while crippling_cond > (mean_t**2 / Area_stringer):
            if web_thickness < data_dictionary["web_maximum_thickness"]:
               web_thickness += data_dictionary["web_minimum_thickness"]*0.01
               Area_stringer = stringer_area(type, web_length, flange_l_prop, web_thickness, \
                                             flange_t_prop)

               if type == "Z" or type == "I" or type == "C":
                   flange_length = web_length * flange_l_prop
                   flange_thickness = web_thickness * flange_t_prop
                   mean_t = web_thickness * (web_length*web_thickness)/Area_stringer + \
                           + flange_thickness * (2*flange_length*flange_thickness)/Area_stringer
               elif type == "T":
                   flange_length = web_length * flange_l_prop
                   flange_thickness = web_thickness * flange_t_prop
                   mean_t = web_thickness * (web_length*web_thickness)/Area_stringer + \
                           + flange_thickness * (flange_length*flange_thickness)/Area_stringer

               load = Mx / (N_stringers*Area_stringer*z_max)
               crippling_cond = ((load*RF) / (B*(E**m * s_yield**(3-m))**(1/3)))**(1/m)

            elif web_thickness >= data_dictionary["web_maximum_thickness"] and \
                 web_length < data_dictionary["web_maximum_length"]:
                 web_length += data_dictionary["web_minimum_length"]*0.01
                 Area_stringer = stringer_area(type, web_length, web_thickness, flange_length, \
                                               flange_thickness)

                 if type == "Z" or type == "I" or type == "C":
                     flange_length = web_length * flange_l_prop
                     flange_thickness = web_thickness * flange_t_prop
                     mean_t = web_thickness * (web_length*web_thickness)/Area_stringer + \
                             + flange_thickness * (2*flange_length*flange_thickness)/Area_stringer
                 elif type == "T":
                     flange_length = web_length * flange_l_prop
                     flange_thickness = web_thickness * flange_t_prop
                     mean_t = web_thickness * (web_length*web_thickness)/Area_stringer + \
                             + flange_thickness * (flange_length*flange_thickness)/Area_stringer

                 load = Mx / (N_stringers*Area_stringer*z_max)
                 crippling_cond = ((load*RF) / (B*(E**m * s_yield**(3-m))**(1/3)))**(1/m)

            elif web_thickness >= data_dictionary["web_maximum_thickness"] and \
                 web_length >= data_dictionary["web_maximum_length"]:
                 print('Maximum stringer dimensions reached while still suffering '+ \
                       'from crippling.')
                 print('Increase maximum allowable dimensions or reduce safety reserve factor.')
                 break

    return web_length, web_thickness, flange_length, flange_thickness
#-------------------------------------------------------------------------------
def skin_length(x, z, boom_index):
    """
    DESCRIPTION: Calculates the skin length to lump its area at the discrete
                 positions that define the wing profile.
    INPUTS:
        * [x]: Python list containing the x coordinates of the section.
        * [z]: Python list containing the y coordinates of the section.
        * [boom_index]: Python list containing the index position where the spar
                        walls are located.
    OUTPUTS:
        * L: Whole skin length of the structure.
    """
    L = 0
    for i in range(len(x)-1):
        L += math.sqrt((x[i+1]-x[i])**2 + (z[i+1]-z[i])**2)

    L += abs(z[boom_index[0]] - z[boom_index[4]]) + \
         abs(z[boom_index[1]] - z[boom_index[3]])

    return L
#-------------------------------------------------------------------------------
def panel_sizing(data_dictionary, y, Y1_list, Vz_list, Mx_list, Vx,
                 Y2_list, My_list):
    """
    DESCRIPTION: Providing the loading conditions in a given wingspan position
                 "y", sizes the skin of the section panels in order to be
                 compliant with the safety factor withstanding the shear
                 loading.
    INPUTS:
        * {data_dictionary}: The dictionary containing all the user input data.
                       All filenames and variables can be obtained from here.
        * y: Span position where loads are desired to be obtained.
        * [Y1_list]: A Python list containing the discrete span positions.
        * [Vz_list]: A Python list containing the vertical shear force at the
                   span positions.
        * [Mx_list]: A Python list containing the bending moment at the span
                   positions.
        * Vx: Drag force.
        * [Y1_list]: A Python list containing the discrete span positions where
                     torsional moment is calculated.
        * [My_list]: A Python list containing the torsional moment at the
                     span positions.
    OUTPUTS: boom_index, Q_star, skin_t, panel_stress
        * [boom_index]: Python list containing the index position where the spar
                        walls are located.
        * [Q_star]: Python list containing the solved circular shear flows.
        * skin_t: Thickness of the skin panels.
        * panel_stress: Maximum stress reached in a panel. This is the stress
                        used for sizing.
    """
    c = chord_variation_law(data_dictionary, y)
    x, z = extract_airfoil(data_dictionary, c)
    Vz, Mx = interpolation_shear_moment(y, Y1_list, Vz_list, Mx_list)
    My = interpolation_torsion(y, Y2_list, My_list)

    x_boom = [x[0]]
    z_boom = [z[0]]
    boom_index = []
    for i in range(1, len(x)):
        if x[i] == 3*c/4 or x[i] == c/5 or x[i] == 0:
            x_boom.append(x[i])
            z_boom.append(z[i])
            boom_index.append(i)

        elif x[i] < 3*c/4 and x[i-1] > 3*c/4:
            x_boom.append(x[i])
            z_boom.append(z[i])
            boom_index.append(i)

        elif x[i] < c/5 and x[i-1] > c/5:
            x_boom.append(x[i-1])
            z_boom.append(z[i-1])
            boom_index.append(i-1)

        elif x[i] > c/5 and x[i-1] < c/5:
            x_boom.append(x[i])
            z_boom.append(z[i])
            boom_index.append(i)

        elif x[i] > 3*c/4 and x[i-1] < 3*c/4:
            x_boom.append(x[i-1])
            z_boom.append(z[i-1])
            boom_index.append(i-1)

    x_boom.append(x[0])
    z_boom.append(z[0]) # First and last points are the same, closed section

    RF = data_dictionary["skin_reserve_factor"]
    skin_t = data_dictionary["skin_minimum_thickness"]

    cond = 0
    while cond < RF:
        # CALCULATION OF THE BOOM AREAS SIMPLIFICATION
        A_boom = skin_length(x, z, boom_index)*(skin_t/1000) / 6 # In meters

        # CALCULATION OF THE MOMENTS OF INERTIA AND STATIC MOMENTS OF INERTIA
        xbg, zbg = gravity_centre(x_boom, z_boom)
        xb_cg, zb_cg = change_to_CG(x_boom, z_boom, xbg, zbg)
        Ixx = 0
        Izz = 0
        Ixz = 0
        for i in range(len(xb_cg)-1):
            Ixx += A_boom * zb_cg[i]**2
            Izz += A_boom * xb_cg[i]**2
            Ixz += A_boom * xb_cg[i]*zb_cg[i]

        Sx = []
        Sz = []
        for i in range(len(xb_cg)-1):
            Sx.append(A_boom * zb_cg[i])
            Sz.append(A_boom * xb_cg[i])

        # CALCULATION OF OPEN SHEAR FLOWS
        q_prime = []
        Sx_cumul = 0
        Sz_cumul = 0
        for i in range(len(xb_cg) - 2):
            Sx_cumul += Sx[i]
            Sz_cumul += Sz[i]
            q_prime.append(-(Vz*Izz - Vx*Ixz)*Sx_cumul / (Ixx*Izz - Ixz**2) \
                           -(Vx*Ixx - Vz*Ixz)*Sz_cumul / (Ixx*Izz - Ixz**2))
        q_prime.append(0)
        q_prime.append(0)
        q_prime.append(0)

        # CALCULATION OF AREAS FOR THE BREATH-BATHO FORMULA
        Omega = []
        for i in range(len(xb_cg)-1):
            Omega.append(abs(xb_cg[i]*zb_cg[i+1] - xb_cg[i+1]*zb_cg[i]) / 2) # In m^2
        Omega.append(abs(xb_cg[1]*zb_cg[5] - xb_cg[5]*zb_cg[1]) / 2)
        Omega.append(abs(xb_cg[2]*zb_cg[4] - xb_cg[4]*zb_cg[2]) / 2)

        # DEFINITION OF THE SYSTEM OF EQUATIONS
        A = np.array([[0., 0., 0.],
                      [0., 0., 0.],
                      [0., 0., 0.]])
        B = np.array([[0.],
                      [0.],
                      [0.]])

        # MOMENT EQUATION
        A[0][0] = -2*Omega[0] - 2*Omega[5] + 2*Omega[6]
        A[0][1] = -2*Omega[1] - 2*Omega[4] - 2*Omega[6] + 2*Omega[7]
        A[0][2] = -2*Omega[2] - 2*Omega[3] - 2*Omega[7]

        B[0][0] = My
        for i in range(5):
            B[0][0] += -2*Omega[i]*q_prime[i]

        # FIRST TWISTING COMPATIBILITY EQUATION
        l = []
        for i in range(len(xb_cg)-1):
            l.append(math.sqrt((xb_cg[i]-xb_cg[i+1])**2 + (zb_cg[i]-zb_cg[i+1])**2))
        l.append(math.sqrt((xb_cg[1]-xb_cg[5])**2 + (zb_cg[1]-zb_cg[5])**2))
        l.append(math.sqrt((xb_cg[2]-xb_cg[4])**2 + (zb_cg[2]-zb_cg[4])**2))

        Area_1 = Omega[0] + Omega[5] - Omega[6]
        Area_2 = Omega[1] + Omega[4] + Omega[6] + Omega[7]
        Area_3 = Omega[2] + Omega[3] - Omega[7]

        A[1][0] = (-l[0]-l[5]-l[6])/Area_1 - l[6]/Area_2
        A[1][1] = l[6]/Area_1 + (l[1] + l[6] + l[4] + l[7])/Area_2
        A[1][2] = -l[7]/Area_2

        B[1][0] = (-q_prime[0]*l[0])/Area_1 + (q_prime[1]*l[1] + q_prime[4]*l[4])/Area_2

        # SECOND COMPATIBILITY EQUATION
        A[2][0] = l[6]/Area_2
        A[2][1] = (-l[1] - l[6] - l[4] - l[7])/Area_2 - l[7]/Area_3
        A[2][2] = l[7]/Area_2 + (l[2] + l[3] + l[7])/Area_3

        B[2][0] = (-q_prime[1]*l[1] - q_prime[4]*l[4])/Area_2 + (q_prime[2]*l[2] + \
                                  q_prime[3]*l[3])/Area_3

        # RESOLUTION OF THE SYSTEM OF EQUATIONS TO OBTAIN THE CIRCULAR CELL SHEAR FLOWS
        Q_star = np.linalg.inv(A).dot(B)

        # CORRECTION OF ALL THE SHEAR FLOWS WITH THE CIRCULAR SHEAR FLOWS
        q = [q_prime[0] - Q_star[0][0]]
        q.append(q_prime[1] - Q_star[1][0])
        q.append(q_prime[2] - Q_star[2][0])
        q.append(q_prime[3] - Q_star[2][0])
        q.append(q_prime[4] - Q_star[1][0])
        q.append(q_prime[5] - Q_star[0][0])
        q.append(q_prime[6] + Q_star[0][0] - Q_star[1][0])
        q.append(q_prime[7] + Q_star[1][0] - Q_star[2][0])

        q_abs = []
        for i in range(len(q)):
            q_abs.append(abs(q[i]))

        max_q = max(q_abs)
        panel_stress = (max_q / (skin_t/1000))*10**(-6)
        cond = data_dictionary["skin_yield_stress"] / panel_stress
        if cond < RF and skin_t < data_dictionary["skin_minimum_thickness"]:
            skin_t += data_dictionary["skin_minimum_thickness"]*0.01
        elif skin_t >= data_dictionary["skin_maximum_thickness"]:
            print("Maximum skin thickness reached while dimensions not compliant" + \
                  "with the safety reserve factor.")
            print("Increase maximum allowable thickness or reduce safety factor.")
            break

    return boom_index, Q_star, skin_t, panel_stress
#-------------------------------------------------------------------------------
def shear_stress_graph_points(data_dictionary, y, Y1_list, Vz_list, Mx_list, Vx,
                              Y2_list, My_list, boom_index, Q_star, skin_t):
    """
    DESCRIPTION: Obtains a list of the shear stress at every panel for a single
                 cell section and for a three-cell sections with the objective
                 of comparing.
    INPUTS:
        * {data_dictionary}: The dictionary containing all the user input data.
                       All filenames and variables can be obtained from here.
        * y: Span position where loads are desired to be obtained.
        * [Y1_list]: A Python list containing the discrete span positions.
        * [Vz_list]: A Python list containing the vertical shear force at the
                   span positions.
        * [Mx_list]: A Python list containing the bending moment at the span
                   positions.
        * Vx: Drag force.
        * [Y1_list]: A Python list containing the discrete span positions where
                     torsional moment is calculated.
        * [My_list]: A Python list containing the torsional moment at the
                     span positions.
        * [boom_index]: Python list containing the index position where the spar
                        walls are located.
        * [Q_star]: Python list containing the solved circular shear flows.
        * skin_t: Thickness of the skin panels.
    OUTPUTS:
        * [q_single_abs]: Python list containing the stress in each panel of the
                          section for a multicell section.
        * [q_multi_abs]: Python list containing the stress in each panel of the
                          section for a single-cell section.
    """
    c = chord_variation_law(data_dictionary, y)
    x, z = extract_airfoil(data_dictionary, c)
    Vz, Mx = interpolation_shear_moment(y, Y1_list, Vz_list, Mx_list)
    My = interpolation_torsion(y, Y2_list, My_list)

    xg, zg = gravity_centre(x, z)
    x_cg, z_cg = change_to_CG(x, z, xg, zg)
    A_lumped = skin_length(x_cg, z_cg, boom_index)*(skin_t/1000) / (len(x)-1)

    # CALCULATION OF THE MOMENTS OF INERTIA
    Ixx = 0
    Izz = 0
    Ixz = 0
    for i in range(len(x_cg)-1):
        Ixx += A_lumped * z_cg[i]**2
        Izz += A_lumped * x_cg[i]**2
        Ixz += A_lumped * x_cg[i]*z_cg[i]

    Sx = []
    Sz = []
    for i in range(len(x_cg)-1):
        Sx.append(A_lumped * z_cg[i])
        Sz.append(A_lumped * x_cg[i])

    # CALCULATION OF OPEN SHEAR FLOWS
    q_prime = []
    Sx_cumul = 0
    Sz_cumul = 0
    for i in range(len(x_cg) - 2):
        Sx_cumul += Sx[i]
        Sz_cumul += Sz[i]
        q_prime.append(-(Vz*Izz - Vx*Ixz)*Sx_cumul / (Ixx*Izz - Ixz**2) \
                       -(Vx*Ixx - Vz*Ixz)*Sz_cumul / (Ixx*Izz - Ixz**2))
    q_prime.append(0)
    q_prime.append(0)
    q_prime.append(0)

    # CALCULATION OF AREAS FOR THE BREATH-BATHO FORMULA
    Omega = []
    for i in range(len(x_cg)-1):
        Omega.append(abs(x_cg[i]*z_cg[i+1] - x_cg[i+1]*z_cg[i]) / 2)

    # CALCULATION OF THE CIRCULAR SHEAR FLOW
    Omega_qp = []
    for i in range(len(x_cg)-1):
        Omega_qp.append(Omega[i] * q_prime[i])

    q_multi_abs = []
    for i in range(boom_index[0]):
        q_multi_abs.append(abs(q_prime[i] - Q_star[0][0]))

    for i in range(boom_index[0], boom_index[1]):
        q_multi_abs.append(abs(q_prime[i] - Q_star[1][0]))

    for i in range(boom_index[1], boom_index[3]):
        q_multi_abs.append(abs(q_prime[i] - Q_star[2][0]))

    for i in range(boom_index[3], boom_index[4]):
        q_multi_abs.append(abs(q_prime[i] - Q_star[1][0]))

    for i in range(boom_index[4], len(q_prime)-2):
        q_multi_abs.append(abs(q_prime[i] - Q_star[0][0]))

    q_single_abs = []
    for i in range(len(q_prime)-2):
        q_single_abs.append(abs(q_prime[i]-Q_star[0][0]))

    for i in range(len(q_multi_abs)):
        q_single_abs[i] = q_single_abs[i]/(skin_t/1000)*10**(-6)
        q_multi_abs[i] = q_multi_abs[i]/(skin_t/1000)*10**(-6)

    return q_single_abs, q_multi_abs
#-------------------------------------------------------------------------------
def skin_buckling_stringer(data_dictionary, y, skin_t, panel_stress):
    """
    DESCRIPTION: Updates the number of stringers taking into account skin
                 buckling calculations with the load at the span section "y".
    INPUTS:
        * {data_dictionary}: The dictionary containing all the user input data.
                       All filenames and variables can be obtained from here.
        * y: Span position where loads are desired to be obtained.
        * skin_t: Thickness of the skin panels.
        * panel_stress: Maximum stress reached in a panel. This is the stress
                        used for sizing.
    OUTPUTS:
        * N_stringers_buckling: The minimum number of stringers needed to avoid
                                skin buckling.
    """
    K_shear = 5.2
    E = data_dictionary["skin_elastic_modulus"]*10**9
    RF = data_dictionary["skin_reserve_factor"]
    b = (skin_t)/1000 * math.sqrt(K_shear*E / (RF*panel_stress*10**6))

    c = chord_variation_law(data_dictionary, y)
    x, z = extract_airfoil(data_dictionary, c)
    L = 0
    for i in range(len(x)-1):
        L += math.sqrt((x[i+1]-x[i])**2 + (z[i+1]-z[i])**2)

    N_stringers_buckling = math.ceil(L/b)

    return N_stringers_buckling
#-------------------------------------------------------------------------------
def neutral_line(data_dictionary, y):
    """
    DESCRIPTION: Calculates the initial and final point of the
    INPUTS:
        * {data_dictionary}: The dictionary containing all the user input data.
                       All filenames and variables can be obtained from here.
        * y: Span position where loads are desired to be obtained.
    OUPUTS:
        * [x_neutral]: Python list containing the x coordinates of the neutral
                       line.
        * [z_neutral]: Python list containing the z coordinates of the neutral
                       line.
    """
    c = chord_variation_law(data_dictionary, y)
    x, z = extract_airfoil(data_dictionary, c)
    xg, zg = gravity_centre(x, z)
    x_cg, z_cg = change_to_CG(x, z, xg, zg)

    Ixx = 0
    Izz = 0
    Ixz = 0
    for i in range(len(x_cg)-1):
        Ixx += z_cg[i]**2
        Izz += x_cg[i]**2
        Ixz += x_cg[i]*z_cg[i]
    # As we are going to divide Ixz by Izz, it is not necessary to multiply by
    # the area, as they are going to be cancelled.

    x_neutral = [min(x_cg), max(x_cg)]
    z_neutral = [Ixz/Izz * x_neutral[0], Ixz/Izz * x_neutral[1]]

    return x_neutral, z_neutral
#-------------------------------------------------------------------------------
def bending_stress_graph(data_dictionary, y, web_l, web_t, skin_t, N_stringers, boom_index, Mx):
    """
    DESCRIPTION:
    INPUTS:
        * {data_dictionary}: The dictionary containing all the user input data.
                       All filenames and variables can be obtained from here.
        * y: Span position where loads are desired to be obtained.
        * web_l: Final web length.
        * web_t: Final web thickness.
        * skin_t: Thickness of the skin panels.
        * N_stringers: Final number of stringers.
        * [boom_index]: Python list containing the index position where the spar
                        walls are located.
        * Mx: The moment at the winspan position "y". It is provided as an
              output in order for its further use and not to be interpolated again.
    OUTPUTS:
        * [s_y]: Python list containing all the direct stresses in order to plot
                 its evolution along the section.
    """
    c = chord_variation_law(data_dictionary, y)
    x, z = extract_airfoil(data_dictionary, c)
    xg, zg = gravity_centre(x, z)
    x_cg, z_cg = change_to_CG(x, z, xg, zg)

    flange_l_prop = data_dictionary["flange_length_proportion"]
    flange_t_prop = data_dictionary["flange_thickness_proportion"]
    type = data_dictionary["stringer_type"]
    A_st = stringer_area(type, web_l, flange_l_prop, web_t, flange_t_prop)/(10*6) #m^2
    A_sk = skin_length(x, z, boom_index)*(skin_t/1000) #m^2

    A_total = N_stringers*A_st + A_sk
    A_lumped = A_total/(len(x)-1)
    Ixx = 0
    Izz = 0
    Ixz = 0
    for i in range(len(x_cg)-1):
        Ixx += A_lumped * z_cg[i]**2
        Izz += A_lumped * x_cg[i]**2
        Ixz += A_lumped * x_cg[i]*z_cg[i]

    s_y = []
    for i in range(len(x_cg)-1):
        s_y.append((Mx*Ixz*x_cg[i]/(Ixx*Izz-Ixz**2) - Mx*Izz*z[i]/(Ixx*Izz-Ixz**2))*10**(-6))

    return s_y
#-------------------------------------------------------------------------------
