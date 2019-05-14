"""
file name: main_script
author: Sacha Gavino
date: March 2019
language: PYTHON 3.7.0
version: 2.1
"""

import os
import subprocess
import sys
import re
from shutil import copytree, copyfile, ignore_patterns, rmtree, move
from tempfile import mkstemp
import numpy as np
import datetime as dt
import fileinput
import pandas as pd
#sys.path.append(os.path.join(os.path.dirname(sys.path[0]), 'PROCEDURE', 'disk_model'))
import disk_structure


#-------------Inputs to change-------------
SOURCE_NAME = "SOURCE_NAME"
NUMBER = "NUMBER"
DATE = "DATE"
Star_mass = "Star_mass"
Ref_radius = "Ref_radius"  
Cut_radius = "Cut_radius"
Inner_radius = "Inner_radius"
Outer_radius = "Outer_radius" 
Tmid_ref = "Tmid_ref"
Tatm_ref = "Tatm_ref" 
Q_exp = "Q_exp"  
Sigma_T = "Sigma_T"  
Nb_points = "Nb_points"  
Surface_density_ref = "Surface_density_ref"   
P_exp = "P_exp"  
NH_to_AV_conversion = "NH_to_AV_conversion"  
Initial_He = "Initial_He"  
Uv_ref = "Uv_ref" 
Small_grains = "Small_grains"  
Big_grains = "Big_grains" 
Dtogas_up = "Dtogas_up" 
Dtogas_in = "Dtogas_in"	
Gamma_exp = "Gamma_exp"	
Transition_altitude = "Transition_altitude"  
Rho_m = "Rho_m"	
D_exp = "D_exp"	
Schmidt_number = "Schmidt_number"	
Alpha = "Alpha"	
S_x	= "S_x"		
Q_c	= "Q_c"		
Settling = "Settling"     	  
A_min = "A_min"   
A_max = "A_max"
Nb_sizes = "Nb_sizes"


#-------------function 1: replace-------------
def sed(pattern, replace, source, dest=None, count=0):
    """Reads a source file and writes the destination file.

    In each line, replaces pattern with replace.

    Args:
        pattern (str): pattern to match (can be re.pattern)
        replace (str): replacement str
        source  (str): input filename
        count (int): number of occurrences to replace
        dest (str):   destination filename, if not given, source will be over written.        
    """
    fin = open(source, 'r')
    num_replaced = count
    if dest:
        fout = open(dest, 'w')
    else:
        fd, name = mkstemp()
        fout = open(name, 'w')
    for line in fin:
        out = re.sub(pattern, replace, line)
        fout.write(out)
        if out != line:
            num_replaced += 1
        if count and num_replaced > count:
            break
    try:
        fout.writelines(fin.readlines())
    except Exception as E:
        raise E
    fin.close()
    fout.close()
    if not dest:
        move(name, source) 

#_________________________________________________#
#                      MAIN                       #
#_________________________________________________#
if __name__ == "__main__":
    with open('../SOURCE/input_parameters/inputs.in') as inputs:
        table = pd.read_table(inputs, sep='    ', index_col=0, comment='#', header=None, lineterminator='\n', engine='python')
    inputs.close()

    today = dt.date.today()
    path_to_source_parameters = "../SOURCE/input_parameters/source_parameters.in"
    path_to_parameters = "../SOURCE/input_parameters/parameters.in"
    path_to_radii = "../SOURCE/input_parameters/radii.in"
    path_to_run_process = "../SOURCE/run_process/run_radius_vanoise.sh"
    #path_to_model = "disk_model/"
    #if os.path.isdir("../RUNS"):
            #rmtree("../RUNS")
    run_directory="RUNS"
    number_of_models = table.shape[1]

    for number in range(1, number_of_models+1, 1):
        name = table.loc['name'][number]
        star_mass = table.loc['star_mass'][number]
        ref_radius = table.loc['ref_radius'][number]
        cut_radius = table.loc['cut_radius'][number]
        inner_radius = table.loc['inner_radius'][number]
        outer_radius = table.loc['outer_radius'][number]
        Tmidplan_ref = table.loc['Tmidplan_ref'][number]
        Tatmos_ref = table.loc['Tatmos_ref'][number]
        q_exp = table.loc['q_exp'][number]
        sigma_T = table.loc['sigma_T'][number]
        nb_points = table.loc['nb_points'][number]
        surface_density_ref = table.loc['surface_density_ref'][number]
        p_exp = table.loc['p_exp'][number]
        nH_to_AV_conversion = table.loc['nH_to_AV_conversion'][number]
        initial_He = table.loc['initial_He'][number]
        uv_ref = table.loc['uv_ref'][number]
        small_grains = table.loc['small_grains'][number]
        big_grains = table.loc['big_grains'][number]
        dtogas_up = table.loc['dtogas_up'][number]
        dtogas_in = table.loc['dtogas_in'][number]
        gamma_exp = table.loc['gamma_exp'][number]
        transition_altitude = table.loc['transition_altitude'][number]
        rho_m = table.loc['rho_m'][number]
        d_exp = table.loc['d_exp'][number]
        schmidt_number = table.loc['schmidt_number'][number]
        alpha = table.loc['alpha'][number]
        s_x = table.loc['s_x'][number]
        q_c = table.loc['q_c'][number]
        settling = table.loc['settling'][number]
        a_min = table.loc['a_min'][number]
        a_max = table.loc['a_max'][number]
        nb_sizes = table.loc['nb_sizes'][number]
        nb_sizes_int = int(float(nb_sizes))
        cutoff = table.loc['cutoff'][number]


        #____change s_param_____#
        sed(SOURCE_NAME, name, path_to_source_parameters, "source_parameters.in")  #"../%s/%s_%s_%s/test_param.in" % (run_directory, today, name, number))
        sed(NUMBER, str(number), "source_parameters.in")
        sed(DATE, str(today), "source_parameters.in")
        sed(Star_mass, star_mass, "source_parameters.in")
        sed(Ref_radius, ref_radius, "source_parameters.in")
        sed(Cut_radius, cut_radius, "source_parameters.in")
        sed(Inner_radius, inner_radius, "source_parameters.in")
        sed(Outer_radius, outer_radius, "source_parameters.in")
        sed(Tmid_ref, Tmidplan_ref, "source_parameters.in")
        sed(Tatm_ref, Tatmos_ref, "source_parameters.in")
        sed(Q_exp, q_exp, "source_parameters.in")
        sed(Sigma_T, sigma_T, "source_parameters.in")
        sed(Nb_points, nb_points, "source_parameters.in")
        sed(Surface_density_ref, surface_density_ref, "source_parameters.in")
        sed(P_exp, p_exp, "source_parameters.in")
        sed(NH_to_AV_conversion, nH_to_AV_conversion, "source_parameters.in")
        sed(Initial_He, initial_He, "source_parameters.in")
        sed(Uv_ref, uv_ref, "source_parameters.in")
        sed(Small_grains, small_grains, "source_parameters.in")
        sed(Big_grains, big_grains, "source_parameters.in")
        sed(Dtogas_up, dtogas_up, "source_parameters.in")
        sed(Dtogas_in, dtogas_in, "source_parameters.in")
        sed(Gamma_exp, gamma_exp, "source_parameters.in")
        sed(Transition_altitude, transition_altitude, "source_parameters.in")
        sed(Rho_m, rho_m, "source_parameters.in")
        sed(D_exp, d_exp, "source_parameters.in")
        sed(Schmidt_number, schmidt_number, "source_parameters.in")
        sed(Alpha, alpha, "source_parameters.in")
        sed(S_x, s_x, "source_parameters.in")
        sed(Q_c, q_c, "source_parameters.in")
        if settling == "yes":
            sed(Settling, "1.000e+00", "source_parameters.in")
        elif settling == "no":
            sed(Settling, "0.000e+00", "source_parameters.in")
        sed(A_min, a_min, "source_parameters.in")
        sed(A_max, a_max, "source_parameters.in")
        sed(Nb_sizes, nb_sizes, "source_parameters.in")
        if cutoff == "yes":
            sed(Cutoff, "1.000e+00", "source_parameters.in")
        elif cutoff == "no":
            sed(Cutoff, "0.000e+00", "source_parameters.in")

        copyfile(path_to_radii, "radii.in")
        copyfile(path_to_parameters, "parameters.in")
        #disk_structure.main()
        #subprocess.call('disk_structure.py')
        os.system('ipython disk_structure.py')
        move("MODEL/", "../%s/%s_%s_%s/" % (run_directory, today, name, number))
        #os.makedirs("../%s/%s_%s_%s" % (run_directory, today, name, number)) # creates the subfolders in RUNS.

        copyfile("source_parameters.in", "../%s/%s_%s_%s/source_parameters.in" % (run_directory, today, name, number))

        for subdir, dirs, files in os.walk("../%s/%s_%s_%s/" % (run_directory, today, name, number)):
            copyfile(path_to_run_process, "%s/run_radius_vanoise.sh" % subdir )
        
        #____write informations about the source_____#
        information=open("../%s/%s_%s_%s/information.txt" % (run_directory, today, name, number), "w+")
        information.write("#---------------------------------------------------------------------------#\n")
        information.write("#--------------------------  SOURCE INFORMATION  ----------------------------\n")
        information.write("#---------------------------------------------------------------------------#\n")
        information.write("name: %s\n" % name)
        information.write("number: %s\n" % number)
        information.write("settling: %s\n" % settling)
        information.write("Mstar (Msun): %s\n" % star_mass)
        information.write("Surface density at ref (g.cm-2): %s\n" % surface_density_ref)
        information.write("Temperature midplan at ref (K): %s\n" % Tmidplan_ref)
        information.write("Temperature atmosphere at ref (K): %s\n" % Tatmos_ref)
        information.write("radius at ref (au): %s\n" % ref_radius)
        information.write("cut radius (au): %s\n" % cut_radius)
        information.write("outer radius (au): %s\n" % outer_radius)
        information.write("p exponent: %s\n" % p_exp)
        information.write("Schmidt number: %s\n" % schmidt_number)
        information.write("UV of reference: %s\n" % uv_ref)
        information.write("#---------------------------  grain information  ----------------------------\n")
        information.write("Viscosity coefficient (or turbulent strength parameter): %s\n" % alpha)
        information.write("S_x: %s\n" % s_x)
        information.write("Q_c: %s\n" % q_c)
        information.write("grain size range (cm): %s -- %s\n" % (a_min, a_max))
        information.write("number of grain sizes: %s\n" % nb_sizes)
        information.write("mass density of the grains: %s\n" % rho_m)
        information.write("#-------------------------------  Comments  ---------------------------------\n")
        
        
        #if os.path.isfile("source_parameters.in"):
            #os.remove("source_parameters.in")
        #if os.path.isfile("radii.in"):
            #os.remove("radii.in")
        #if os.path.isfile("parameters.in"):
            #os.remove("parameters.in")



    