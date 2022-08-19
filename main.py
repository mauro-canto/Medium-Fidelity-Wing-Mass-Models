# IMPORT ALL CREATED FUNCTIONS
from functions import *
import matplotlib.pyplot as plt
import numpy as np

def run_calculation():

    # READING OF THE INPUT DATA:
    data_dictionary = read_data_input("Data_Input.txt")
    extract_chords(data_dictionary)

    # EXTRACTION OF LOADS FROM AVL FILES
    Y1_list, Vz_list, Mx_list = extract_shear_moment(data_dictionary)
    Vx = extract_drag(data_dictionary)
    Y2_list, My_list = extract_torsion(data_dictionary)

    y = max(min(Y1_list), min(Y2_list))
    n_rib = 1
    y_list = [y]
    N_list = []
    A_list = []

    while y < max(max(Y1_list), max(Y2_list)):
        print('                                         -----------------------')
        print('                                        |                       |')
        print('                                        | CALCULATION OF RIB ' + str(n_rib) + ' |')
        print('                                        |                       |')
        print('                                         -----------------------')
        print('\n')
        # STRINGER PRELIMINARY SIZING
        print(' -----------------')
        print('| STRINGER SIZING |')
        print(' -----------------')
        print('Stringer optimization for obtention of preliminary dimensions:')
        print('-----------------------------------')
        Mx, z_max, N_stringers, web_l, web_t, flange_l, flange_t = \
        stringer_preliminary_sizing(data_dictionary, y, Y1_list, Vz_list, Mx_list)
        print('Number of stringers: ' + str(N_stringers))
        print('Web length: ' + str(round(web_l, 2)) + ' mm')
        print('Web thickness: ' + str(round(web_t, 2)) + ' mm')
        print('Flange length: ' + str(round(flange_l, 2)) + ' mm')
        print('Flange thickness: ' + str(round(flange_t, 2)) + ' mm')
        print('\n')

        # STRINGER CRIPPLING
        print('Stringer crippling check:')
        print('-----------------------------------')
        web_length, web_thickness, flange_length, flange_thickness = \
        stringer_crippling(data_dictionary, Mx, z_max, N_stringers, web_l, web_t)
        print('\n')
        print('Final stringer dimensions:')
        print('-----------------------------------')
        print('Web length: ' + str(round(web_length, 2)) + ' mm')
        print('Web thickness: ' + str(round(web_thickness, 2)) + ' mm')
        print('Flange length: ' + str(round(flange_length, 2)) + ' mm')
        print('Flange thickness: ' + str(round(flange_thickness, 2)) + ' mm')
        print('\n')

        # PANEL SIZING
        print(' -------------')
        print('| SKIN SIZING |')
        print(' -------------')
        boom_index, Q_star, skin_t, panel_stress = panel_sizing(data_dictionary, y, \
        Y1_list, Vz_list, Mx_list, Vx, Y2_list, My_list)

        q_single_abs, q_multi_abs = shear_stress_graph_points(data_dictionary, y, \
        Y1_list, Vz_list, Mx_list, Vx, Y2_list, My_list, boom_index, Q_star, skin_t)

        # SKIN BUCKLING
        N_stringers_buckling = skin_buckling_stringer(data_dictionary, y, skin_t, panel_stress)

        print('Skin panel thickness: ' + str(skin_t) + ' mm')
        print('\n')
        print('Skin buckling check:')
        print('-----------------------------------')
        print('Number of stringers to avoid skin buckling: ' + str(max(N_stringers, N_stringers_buckling)))
        print('\n')

        s_y = bending_stress_graph(data_dictionary, y, web_length, web_thickness, skin_t, \
                                   max(N_stringers, N_stringers_buckling), boom_index, Mx)

        # GRAPHS
        print(' ----------------')
        print('| GRAPH DIAGRAMS |')
        print(' ----------------')
        # WING SECTION
        c = chord_variation_law(data_dictionary, y)
        x, z = extract_airfoil(data_dictionary, c)
        xg, zg = gravity_centre(x, z)
        x_cg, z_cg = change_to_CG(x, z, xg, zg)
        plt.figure(figsize = (15,4))
        plt.plot(x_cg, z_cg, 'b', label = "Wing section")
        plt.plot([x_cg[boom_index[0]], x_cg[boom_index[4]]], \
                 [z_cg[boom_index[0]], z_cg[boom_index[4]]], 'b')
        plt.plot([x_cg[boom_index[1]], x_cg[boom_index[3]]], \
                 [z_cg[boom_index[1]], z_cg[boom_index[3]]], 'b')

        plt.plot(0, 0, 'ro', label = "Centre of gravity")

        x_neutral, z_neutral = neutral_line(data_dictionary, y)
        plt.plot(x_neutral, z_neutral, 'r--', label = 'Neutral line')

        plt.xlabel("X")
        plt.ylabel("Z")
        plt.title("Wing section")
        plt.legend()
        plt.show()

        # BENDING AND SHEAR
        fig, axs = plt.subplots(1,2,figsize = (15,4))

        plt.subplot(1,2,1)
        x1 = np.linspace(1, len(x)-1, len(x)-1)
        plt.plot(x1, s_y, 'g')
        plt.title('Direct stress distribution')
        plt.xlabel('Coordinate along wing surface')
        plt.ylabel('$\sigma_y [MPa]$')

        plt.subplot(1,2,2)
        plt.plot(x1, q_multi_abs, 'g')
        plt.title('Shear stress distribution')
        plt.xlabel('Coordinate along wing surface')
        plt.ylabel('$τ_{xz} [MPa]$')

        plt.show()

        # RIB PLACING
        chord_root = data_dictionary["chord_root"]
        chord_tip = data_dictionary["chord_tip"]
        c = chord_variation_law(data_dictionary, y)
        y_min = min(Y1_list)
        y_max = max(Y2_list)
        plt.figure(figsize = (15,4))
        plt.plot([y_min, y_min], [0, chord_root],'b')
        plt.plot([y_min, y_max], [0, 0],'b')
        plt.plot([y_max, y_max], [0, chord_tip],'b')
        plt.plot([y_min, y_max], [chord_root, chord_tip],'b')
        plt.plot([y, y], [0, c], 'r')
        plt.xlabel('Span position (Y)')
        plt.ylabel('X')
        plt.title('Rib position')
        plt.show()

        # CALCULATION OF NEXT RIB DUE TO BUCKLING
        flange_l_prop = data_dictionary["flange_length_proportion"]
        flange_t_prop = data_dictionary["flange_thickness_proportion"]
        type = data_dictionary["stringer_type"]
        RF = data_dictionary["skin_reserve_factor"]
        E = data_dictionary["skin_elastic_modulus"]*10**9
        A_s = stringer_area(type, web_length, flange_l_prop, web_thickness, \
                                      flange_t_prop)

        N_list.append(max(N_stringers, N_stringers_buckling))
        A_list.append(A_s)

        K_compression = 3.6
        K_shear = 5.2
        F_compression = Mx / (N_stringers*A_s*z_max*1000) # Pa
        F_shear = panel_stress*10**6 # Pa

        l_shear = (skin_t/1000) * math.sqrt(K_shear*E/(RF*F_shear))
        l_compression = (skin_t/1000) * math.sqrt(K_compression*E/(RF*F_compression))

        l = min(l_shear, l_compression)

        y += l
        y_list.append(y)
        n_rib += 1
        print('\n')
        print('\n')

    print('                                         -----------------------')
    print('                                        |                       |')
    print('                                        | WHOLE RIB ARRANGEMENT |')
    print('                                        |                       |')
    print('                                         -----------------------')
    plt.figure(figsize = (15,4))
    plt.plot([y_min, y_min], [0, chord_root],'b')
    plt.plot([y_min, y_max], [0, 0],'b')
    plt.plot([y_max, y_max], [0, chord_tip],'b')
    plt.plot([y_min, y_max], [chord_root, chord_tip],'b')
    plt.xlabel('Span position (Y)')
    plt.ylabel('X')
    plt.title('Rib position')
    del(y_list[-1])
    for i in range(len(y_list)):
        c = chord_variation_law(data_dictionary, y_list[i])
        plt.plot([y_list[i], y_list[i]], [0, c], 'r')
    plt.xlabel('Span position (Y)')
    plt.ylabel('X')
    plt.title('Rib position')
    plt.show()

    # MASS ESTIMATION
    c_root = data_dictionary["chord_root"]
    c_tip = data_dictionary["chord_tip"]
    x_root, z_root = extract_airfoil(data_dictionary, c_root)
    x_tip, z_tip = extract_airfoil(data_dictionary, c_tip)
    L_root = skin_length(x_root, z_root, boom_index)
    L_tip = skin_length(x_tip, z_tip, boom_index)
    rho_skin = data_dictionary["skin_material_density"]
    rho_stringer = data_dictionary["stringer_material_density"]

    b = data_dictionary["span"]
    d = math.sqrt((b/2)**2 + (3/4*c_root - 3/4*c_tip)**2)

    mass = (skin_t/1000)*(L_root+L_tip)*b*rho_skin + max(N_list)*(max(A_list)*10**(-6))*2*d*rho_stringer
    print('\n')
    print('\n')
    print('                                      ----------------------------')
    print('                                     |                            |')
    print('                                     | STRUCTURAL MASS ESTIMATION |')
    print('                                     |                            |')
    print('                                      ----------------------------')
    print('Total structural mass: ' + str(round(mass,2)) + ' kg')
