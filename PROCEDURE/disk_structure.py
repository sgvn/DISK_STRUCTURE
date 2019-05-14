#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
file name: disk_structure
author: Sacha Gavino
date: March 2019
language: PYTHON 3.7
__________________________________________________________________________________________
short description:  The main goal of this code is to create the output tables
                    1D_grain_sizes_R.in that has to be used by Nautilus. 
                    The model used herein is the model described in the documentation.
                    Please refer to README.md for a more detailed description.
__________________________________________________________________________________________

INPUT:

 constants.py: gives a banch of physical and chemical constants usefull for the code.

 radius.in: each line gives the radius in au that will be computed.

 source_parameters.in: gives the input physical parameters of the source.

 
OUTPUT:

 1D_grain_sizes.in: gives grain information at each spatial point. Each line is a 
 spatial point. There is one file per radius. We store them in output folder.
 
 1D_static.dat

 parameters.in

 information.in: gives information on the models computed.  

 Plots... 
 
 
IMPORTANT: this code can't be executed if at least one of the input files is missing. 
__________________________________________________________________________________________
"""


import os
import glob
import sys
import re
import fileinput
from shutil import copytree, ignore_patterns, rmtree, move, copyfile, copy
from tempfile import mkstemp

import numpy as np
import math
import argparse

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma']

import constants as cst

#-------------radius-------------
R01 = "R01"
R02 = "R02"
R03 = "R03"
R04 = "R04"
R05 = "R05"
R06 = "R06"
R07 = "R07"
R08 = "R08"
R09 = "R09"
R10 = "R10"
R11 = "R11"
R12 = "R12"
R13 = "R13"
R14 = "R14"
R15 = "R15"
R16 = "R16"
R17 = "R17"

#-------------spatial resolution-------------
nb_points = "nb_points"

#-------------disk parameters-------------
star_mass = "star_mass"
dust_mass = "dust_mass"
surface_density_ref ="surface_density_ref"
Tmidplan_ref = "Tmidplan_ref"
Tatmos_ref = "Tatmos_ref"
ref_radius = "ref_radius"
cut_radius = "cut_radius"
inner_radius = "inner_radius"
outer_radius = "outer_radius"
sigma_T = "sigma_T"
q_exp = "q_exp"
rho_m = "rho_m"
dtogas_in = "dtogas_in"
dtogas_out = "dtogas_out"
d_exp = "d_exp"
p_exp = "p_exp"
gamma_exp = "gamma_exp"
schmidt_number = "schmidt_number"	
alpha = "alpha"
s_x = "s_x"
q_c = "q_c"
settling = "settling"
dtogas_up = "dtogas_up"
initial_He = "initial_He"
transition_altitude = "transition_altitude"
uv_ref = "uv_ref"
nH_to_AV_conversion = "nH_to_AV_conversion"	
UV_FLUX = "UV_FLUX"

#-------------sizes parameters-------------
a_min = "a_min"
a_max = "a_max"
nb_sizes = "nb_sizes"
small_grains = "small_grains"
big_grains = "big_grains"
NBR_GRAINS = "NBR_GRAINS"

#_______________________________________________#
#                   FUNCTION                    #
#_______________________________________________#
#-------------function 1: get source parameter values-------------
def source_value(param):
    source=open("source_parameters.in","r") #extract values from source_parameters.in
    for line in source:
        line = line.partition('#')[0]
        if param in line:
            parts = line.split()
    source.close()
    return float(parts[1])
	
#-------------function 2: get radii values-------------
def rad_value(rad):
    radii=open("radii.in","r") #read values from radius.in
    for line in radii:
        if rad in line:
            if line.strip():
                parts = line.split()
    radii.close()
    return float(parts[1])
	
#-------------function 3: get number of grain species-------------
#def nb_sizes():
	#non_blank_count = 0
	#size=open("grainsizes.in","r") #read lines of grainsizes.in
	#for line in size:
		#if line.strip():
			#non_blank_count += 1
	#size.close()
	#return non_blank_count
	
	
#-------------function 4: plot profile-------------
def plot(z, H, size, radius, array1, array2):
	fig, ax = plt.subplots(figsize=(12, 10))
	plt.title('density profile: size: %scm and rad: %sau' % (size, radius))
	plt.xlabel('z', fontsize = 17)
	plt.ylabel('density (g/cm-3)', fontsize = 17)
	
	plt.xlim(0, 4)
	plt.ylim(1e-23, 1e-10)		
	ax.set_yscale('log')
	ax.plot(z/H, array1, label='gas')
	ax.plot(z/H, array2, label='dust') #,'.', markersize=4
	
	ax.grid()
	ax.legend()
	plt.show()
	
#-------------function 5: replace-------------
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


#_______________________________________________#
#                     CLASS                     #
#_______________________________________________#
#-------------Class: gas parameters-------------	
class Gas:
    """
	This Class gives the physical quantities related to the gas by using the model.
	Basically, the gas is given by Hydrogren nuclei. Equation references are related to the documentation.
	Instance methods:
		A) initializer (scale height at reference radius)
        B) dens_mid(self, sigma, Hg)
		C) n_H(self, H, r, z)
        D) n_H2(self, H, r, z)
		E) H(self, q_exp, radius, R_ref)
        F) Hatmos_ref(self, q_exp, radius, R_ref)
		G) temperature(self, Tmid_ref, radius, R_ref, q_exp)
        H) temp_atmos(self, Tmid_ref, radius, R_ref, q_exp)
        I) vertical_temperature(self, Tmidplan, Tatmos, z, z0, sigma)
		J) surdens(self, radius) 
		K) surdens_big(self, radius)
        L) omega2(self, radius, solar_mass)
		M) H_density(self, rad) --> works only with Nautilus files
		N) altitude(self, rad) --> works only with Nautilus files
		O) altitude2(slef, rad)
        P) dz(self, z)
        Q) dtogas_mid(self, nb_points, transition_altitude, dtogas_in, dtogas_up, z, nH2)
        R) he_ab(self, dtogas_up, He_abundance)
        S) gtodm(self, nb_points, H, nH, transition_altitude, big_grains, small_grains, dtogas_up, nH_to_AV_conversion, mass_cons, z)
        T) UV_factor(self, UV_ref, ref_radius, radius)
    """

	
	
    """ A)
	Equation (6). Initializer. Gives the gas scale height. Unit: cm.
	INIT function for the class Gas. 
	"""
    def __init__(self, Tmid_ref, R_ref, star_mass):
        H_ref = np.sqrt((cst.kb*Tmid_ref*(R_ref*cst.autocm)**3)/(cst.mu*cst.amu*cst.G*\
        star_mass*cst.Msun))
        self.H_ref = H_ref


    """ B)
    Equation (2). Gas density in the midplane. Unit: cm-3.
    Arguments:
        -sigma: gas surface density.
        -Hg: gas scale height.
    """
    def dens_mid(self, sigma, Hg):
        return (sigma)/(Hg*np.sqrt(2*np.pi)*cst.mu*cst.amu)


    """ C)
    Equation (2). Gas number density. Unit: part/cm^3. ISOTHERMAL ONLY.
    WARNING: We divide here by mass to have the number density instead of mass density.
    Arguments: 
	    -H: gas scale height.
	    -r: radii at which the gas number density will be computed.
	    -z: altitude at which the gas number density will be computed.
    """		
    def n_H(self, surf_dens_ref, H, radius, R_ref, p_exp, z):#doesn't give the same values\
	    #as in Nautilus. We create H_density to extract values from Nautilus model.
	    return 2.*(surf_dens_ref/(cst.mu*cst.amu*H*np.sqrt(2.*np.pi)))*(radius/R_ref)**\
	    (-p_exp)*np.exp(-(z)**2./(2.*H**2.))


    """ D)
    Equation (2). Gas number density. Unit: part/cm^3.
    WARNING: We divide here by mass to have the number density instead of mass density.
    Arguments: 
        -H: gas scale height.
        -r: radii at which the gas number density will be computed.
        -z: altitude at which the gas number density will be computed.
    """		
    def n_H2(self, dens_mid, T, nb_points, omega2, z): #doesn't give the same values
        log = 0
        rho = 0
        log_list = []
        log_list.append(0)
        rho_scaled = 0
        for i in range(1, nb_points, 1):
            log = log_list[i-1] - (np.log(T[i]) - np.log(T[i-1])) - ((cst.mu*\
            cst.amu*omega2*z[i]*(z[i] - z[i-1]))/(cst.kb*T[i]))
            log_list.append(log)
        log_rho = np.asarray(log_list)
        rho = np.exp(log_rho)
        
        maximum=np.amax(rho)
        rho_scaled = (rho*dens_mid)/maximum
        return rho_scaled


    """ E)
    Equation (4). Gas scale height using the power law. Unit: cm.
    Argument: 
    	-radius: radii at which H will be computed.
    """	
    def H(self, q_exp, radius, R_ref):
	    h = (3./2.) - (q_exp/2.)
	    return self.H_ref*(radius/R_ref)**h


    """ F)
	Equation not yet written. Temperatue in the atm. Unit: Kelvin   
    Argument:
		-radius: radii at which the temperature will be computed.
        -Tatm_ref: temperature of reference
	"""
    def Hatmos_ref(self, Tatmos_ref, R_ref, star_mass):
        return np.sqrt((cst.kb*Tatmos_ref*(R_ref*cst.autocm)**3)/(cst.mu*cst.amu*cst.G*\
        star_mass*cst.Msun))
	

    """ G)
    Equation (5). Temperatue in the midplane. Unit: Kelvin
    Argument: 
	    -radius: radii at which the temperature will be computed. 
        -Tmid_ref: temperature of reference
    """
    def temperature(self, Tmidplan_ref, radius, R_ref, q_exp):
        return Tmidplan_ref*(radius/R_ref)**(-q_exp)
		
		
    """ H)
	Equation not yet written. Temperatue in the atm. Unit: Kelvin   
    Argument:
		-radius: radii at which the temperature will be computed.
        -Tatmos_ref: temperature of reference in atmosphere
	"""
    def temp_atmos(self, Tatmos_ref, radius, R_ref, q_exp):
        return Tatmos_ref*(radius/R_ref)**(-q_exp)


    """ I)
    Equation not yet written. Temperature vertical profile. Unit: Kelvin.
    Using the definition by Williams and Best (2014)
    Arguments:
        -Tmid_ref: temperature of reference in midplane
        -Tatmos_ref: temperature of reference in atmosphere
        -z: elevation avove the midplane
    """
    def vertical_temperature(self, Tmidplan, Tatmos, z, z0, sigma):
        return Tmidplan+(Tatmos-Tmidplan)*np.sin((np.pi*z)/(2*z0))**(2*sigma)
        

    """ J)
	Equation (7). Surface density of the gas. Unit: g/cm-2
	Argument: 
		-radius: radii at which the surface density will be computed.
    """
    def surdens(self, surf_dens_ref, radius, R_ref, p_exp):
        return 2*surf_dens_ref*(radius/R_ref)**(-p_exp)  #2 or not??
		
		
    """ K)
    Equation (7b). Surface density of the gas with exponential edge. Unit: g/cm-2
    Argument: 
	    -radius: radii at which the surface density will be computed.
    """	
    def surdens_big(self, surf_dens_ref, radius, R_ref, gamma_exp, R_cut): 
	    return 2*surf_dens_ref*(radius/R_ref)**(-gamma_exp)* np.exp(-(radius/R_cut)**\
	    (2-gamma_exp))


    """ L)
    Equation (8). Keplerian angular velocity.
    Argument:
        -radius: radii at which the velocity is calculated
        -solar_mass: mass of central object.
    """
    def omega2(self, radius, star_mass):
        return (cst.G*star_mass*cst.Msun)/(cst.autocm*radius)**3


    """ M)
    Gives the gas number density. Like method F) but instead the method here extracts 
    the values of the gas number density from Nautilus input file and store them in a list
    densities_list. By definition the list size is the number of spatial points (altitudes
    z computed). We created this method because the returns of method F) don't match the 
    Nautilus values. This method is meant to be temporary.
    Argument: 
	    -rad: radii at which the gas number density is computed. WARNING!! "rad" should be 
	    written like in the ".txt" files from which we extract values i.e. 040, 080 or 200 
	    for examples.
    """	
    def H_density(self, rad):#!"rad" should be written like. 040, 080 or 200 for examples!
	    densities_list = []
	    densities=open("../physical_structure/%sAU.txt" % rad,"r")
	    for line in densities:
		    parts = line.split()
		    densities_list.append(parts[1])
	    densities.close()
	    dens = map(float, densities_list[12:])
	    return dens
		
		
    """ N)
    Gives the altitudes at each spatial point computed by Nautilus. This method extracts 
    the values of the altitude from Nautilus input file and store them in a list z_list.
	By definition the list size is the number of spatial points (altitudes z computed).
	This can be executed only side by side with Nautilus. Use method I) for else.
	Argument: 
		-rad: radii at which the altitudes are computed. WARNING!! "rad" should be written
		 like in the ".txt" files from which we extract values i.e. 040, 080 or 200 for 
		 examples.
	"""	
    def altitude(self, rad):#!! "rad" should be written like. 040, 080 or 200 for examples!			 
        z_list = []
        zs=open("../physical_structure/%sAU.txt" % rad,"r")
        for line in zs:
            parts = line.split()
            z_list.append(parts[0])
        zs.close()
        z = map(float, z_list[12:])
        return z # z in au!
		
		
    """ O)
	Gives the altitudes at each spatial point. This method stores the altitude in a list.
	By definition the list size is the number of spatial points (altitudes z computed).
	Argument: 
		-nb_points: number of altitudes computed.
		-H_g: scale height of the gas at the considered radius.
	"""	
    def altitude2(self, nb_points, H_g):		 
        z_list = []
        z = 0
        for i in range(0, nb_points, 1):
            z = (1. - (2.*i/(2.*nb_points - 1.)))*4.*H_g
            z_list.append(z)
        z_list_np = np.asarray(z_list)
        return z_list_np


    """ P)
    Gives the difference in altitude between two consecutive spatial points.
    usefull to compute the Av.
    Arguments:
        -z: elevation.
        -nb_points: number of spatial points computed in Nautilus.
    """
    def diff(self, z, nb_points):
        dz = [0]
        z_b = np.append(z, 0.)
        for spatial_pt in range(0, nb_points-1,1):
            diff = (z_b[spatial_pt] - z_b[spatial_pt+1])
            dz.append(diff)
        return dz

    """ Q)
    Compute the mass conservation to compute the correct gtodn as a function
    of elevation. from Wakelam and Majumdar (2017)
    Arguments:
        -transition_altitude: elevation from where small grains become big.
        -z: elevation
    """
    def dtogas_mid(self, nb_points, transition_altitude, dtogas_in, dtogas_up, z, nH2):
        mass1 = mass2 = 0
        for i in range(1, nb_points, 1):
            if (z[i]/H) > transition_altitude:
                mass1 = mass1 + nH2[i]*cst.mu*cst.amu
            elif (z[i]/H) <= transition_altitude:
                mass2 = mass2 + nH2[i]*cst.mu*cst.amu
        return (dtogas_in*(mass1+mass2) - dtogas_up*mass1)/mass2


    """ R)
    disk_structure.f90. Recompute dust to gas massratio to remove He. In the following, 
    dtogm is used as a H/dust mass ratio.
    from Wakelam and Majumdar (2017).
    Argument:
	    -He_ab: Ab of Helium used in the model.
    """
    def he_ab(self, dtogas_up, He_abundance):
        return dtogas_up*(1. + 4.*He_abundance)



    """ S)
    disk_structure.f90. Return four arrays: grain_size, 1/abundance, AV/NH, Av.
    Used in the case of one size distribution only. The Av is defined here. Notice
    that the class Grains also has a definition of the Av but only in the case of
    vertically isothermal disks for the moment, so it is not used to compute the 
    1D_static.dat for the moment.
    Arguments:
        -transition_altitude: elevation from which small grains become big.
        -grain_r: grain size.
        -rho_m: material grain density.
        -dtogas_up: dtogas mass ratio above transition.
        -NH_to_AV_conversion: conversion factor.
        -mass_cons: correct dtogas (see mass_cons method in Grains class).
    """
    def gtodm(self, nb_points, H, nH, transition_altitude, rho_m, big_grains, small_grains, dtogas_up, nH_to_AV_conversion, dtogas_mid, z):
        grain_r = []
        gtodm = []
        AV_NH = []
        Av_z = []
        grain_r.append(small_grains)
        gtodm.append((4.*np.pi*rho_m*small_grains**3)/(3.*dtogas_up*cst.amu))
        AV_NH.append((1/nH_to_AV_conversion)*(dtogas_up/1e-2)*(1e-5/small_grains))
        Av_z.append(0)
        
        for i in range(1, nb_points, 1):
            if (z[i]/H) > transition_altitude:
                grain_r.append(small_grains)
                gtodm.append((4.*np.pi*rho_m*small_grains**3)/(3.*dtogas_up*cst.amu))  
                AV_NH.append((1/nH_to_AV_conversion)*(dtogas_up/1e-2)*(1e-5/small_grains))
            elif (z[i]/H) <= transition_altitude:
                grain_r.append(big_grains)
                gtodm.append((4.*np.pi*rho_m*big_grains**3)/(3.*dtogas_mid*cst.amu))
                AV_NH.append((1/nH_to_AV_conversion)*(dtogas_mid/1e-2)*(1e-5/big_grains))
            Avz = Av_z[i-1] + nH[i]*AV_NH[i]*(z[i-1]-z[i])*(1e-5/grain_r[i])
            Av_z.append(Avz)
        grain_rnp = np.asarray(grain_r)
        gtodm_np = np.asarray(gtodm)
        AV_NH_np = np.asarray(AV_NH)
        Av_z_np = np.asarray(Av_z)
        return grain_rnp, gtodm_np, AV_NH_np, Av_z_np
                

    """ T)
    disk_structure.f90. Computation of the UV irradiation factor at the requested radius.
    the UV is divided by two because since we do not make a full 3D transfer we assume 
    that only half of the photons irradiated from the star diffuse into the disk, the 
    other half diffuse in the opposite direction. The UV flux coming from the star is 
    assumed to have the same spectrum as the interstellar DRAINE field multiplied by a 
    certain amount and diluted as 1/r^2.
    Arguments:
        -transition_altitude: elevation from which small grains become big.
    """
    def uv_factor(self, UV_ref, ref_radius, radius, Hg):
        return UV_ref/(2*((radius/ref_radius)**2 + ((4*Hg)/(ref_radius*cst.autocm))**2))


	
#-------------Class: grain parameters-------------	
class Grains:
    """
    This Class gives the physical quantities related to the grains using the model.
    Equation references are related to the documentation.
    Instance methods:
	    A) initializer (store the max and min value of grain size ranges.)
	    B) mass_density(self, d_exp, rho_m, C_cst)
	    C) C(self)
	    D) single_surf_dens(self, single_surdens0, radius, R_ref, p_exp)
	    E) single_surf_dens_big(self, single_surdens0, radius, R_ref, gamma_exp,cut_radius)
	    F) rho_d(self, sigma, H_d, z)
	    G) H_d(self, H_g, T_s0, Sc, alpha)
	    H) stopping_time(self, a_av, rho_g, H_g)
	    I) stopping_time0(self, a_av, sigma_g)
	    J) Ts_t(self, alpha, Sc, s_x)
	    K) a_t(self, Ts_t, sigma_ref, rho_m, Rc, ref_radius, p_exp)
	    L) single_surf_dens0(self, dtogas, fraction, gas_surf_dens0)
	    M) average(self, d_exp)
        N) single_mass(self, rho_m, size)
        O) q_ext(self, A, lambd, size_av, q_c)
        P) Av(self, n_d, Q_ext, size_av, z)
    """

    

    """ A)
    Store the grain size maximum values of the ith interval into variables.
    Initializer for the class Grains. 
    """
    """
    def __init__(self, i): # we use the i index to loop on the intervals in grainsizes.in.
        grainsizes=open("grainsizes.in","r")
        lines=grainsizes.readlines()

        parts = lines[i].split()
        min = parts[1]
        max = parts[2]
        self.minf = float(min)
        self.maxf = float(max)
        grainsizes.close()
    """

    def __init__(self, i, a_min, a_max, nb_sizes):
        a = np.logspace(np.log10(a_min), np.log10(a_max), nb_sizes+1)
        self.minf = float(a[i])
        self.maxf = float(a[i+1])
	
	
    """ B)
	Equation (17). mass of dust in the ith interval. Unit: g
	Argument: 
		-d_exp: size distribution exponent.
		-rho_m: material density of the grains.
		-C_cst: normalization constant in the MRN equation. !!! This must be an output !!!
	"""
    def mass_density(self, d_exp, rho_m, C_cst):#this is divided by the gas number density.
	    return ((4.*np.pi)/3.*(4.-d_exp))*rho_m*C_cst*(self.maxf**(4.-d_exp) - \
	    self.minf**(4.-d_exp))
	
	
    """ C)
	Equation (19). normalization constant in the MRN equation
	Argument:
		-single_surdens0: surface density at reference radius.
		-radius: radii at which the surface density will be computed.
		-R_ref: reference radius.
		-p_exp: surface density exponent. 
	"""	
    def C(self, dtogas, rho_m, d_exp, amax, amin):
	    return (3./(4.*np.pi))*((dtogas*cst.mu*cst.proton_mass_g)/rho_m)*(4 - d_exp)*(1/(amax**(4-d_exp) - amin**(4-d_exp)))
		
		
    """ D)
	Equation (20). Surface density of single grain species. Unit: g/cm^2
	Argument:
		-single_surdens0: surface density at reference radius.
		-radius: radii at which the surface density will be computed.
		-R_ref: reference radius.
		-p_exp: surface density exponent. 
	"""	
    def single_surf_dens(self, single_surdens0, radius, R_ref, p_exp):
        return single_surdens0*(radius/R_ref)**(-p_exp)
		
		
    """ E)
	Equation (20b). Surface density with exponential edge. Unit: g/cm^2
	Argument:
		-single_surdens0: surface density at reference radius.
		-radius: radii at which the surface density will be computed.
		-R_ref: reference radius.
		-gamma_exp: surface density exponent.
		-cut_radius: tapered radius. 
	"""	
    def single_surf_dens_big(self, single_surdens0, radius, R_ref, gamma_exp, cut_radius):
	    return single_surdens0*(radius/R_ref)**(-gamma_exp)*np.exp(-(radius/cut_radius)\
	    **(2-gamma_exp))
	
	
    """ F)
	Equation (23). mass density of single grain species. Unit: g/cm^3
	Argument: 
		-sigma: surface density of single grain species.
		-H_d: dust scale height of single grain species.
		-z: altitudes above midplane.
	"""	
    def rho_d(self, single_sigma_d, H_d, z):
	    return (single_sigma_d/(np.sqrt(2*np.pi)*H_d))*np.exp(-(z**2)/(2*H_d**2))
				
	
    """ G)
	Equation (26). dust scale height of single grain species. Unit: cm
	Argument: 
		-H_g: gas scale height.
		-T_so: dimensionless stopping time in the midplane.
		-Sc: Schmidt number.
		-alpha: viscosity coefficient.
	"""	
    def H_d(self, H_g, T_s0, Sc, alpha):
        return H_g/(np.sqrt(1 + T_s0*(Sc/alpha)))
	
	
    """ H)
	Equation (30). Dimensionless stopping time.
	Argument: 
		-a_av: averaged grain size of the ith interval.
		-rho_m: material density of the grains.
		-rho_g: mass density of the gas.
		-H_g: gas scale height at considered radius.
	"""	
    def stopping_time(self, a_av, rho_m, rho_g, H_g):
	    return (a_av*rho_m)/(rho_g*H_g)
		
		
    """ I)
	Equation (31). Dimensionless stopping time in the midplane.
	Argument: 
		-a_av: averaged grain size of the ith interval.
		-rho_m: material density of the grains.
		-sigma_g: gas surface density.  
	"""
    def stopping_time0(self, a_av, rho_m, sigma_g):
        return (np.sqrt(2*np.pi)*a_av*rho_m)/sigma_g
	
	
    """ J)
	Equation (33). Dimensionless stopping time in midplane at which the settling factor is
				   x % smaller than the asymptotic value for small grains. 
	Argument: 
		-alpha: viscosity coefficient.
		-Sc: Schmidt number.
		-s_x: settling factor between big grains and small grains.
	"""
    def Ts_t(self, alpha, Sc, s_x):
        return (alpha/Sc)*((1. - s_x**2.)/s_x**2.)
		
		
    """ K)
	Equation (34). Transition size between small and big grains.
	Argument: 
		-Ts_t: transition dimensionless stopping time in midplane.
		-sigme_ref: gas surface density at reference radius.
		-rho_m: material mass density of the grains.
		-Rc: cut radius at which you consider the transition size.
		-ref_radius: reference radius.
		-p_exp: surface density exponent.
	"""
    def a_t(self, Ts_t, sigma_ref, rho_m, cut_radius, ref_radius, p_exp):
	    return ((Ts_t*sigma_ref)/(np.sqrt(2*np.pi)*rho_m))*(cut_radius/ref_radius)**(-p_exp)
	
	
    """ L)
	Equation (43). Surface density of single sized grains at reference radius. Unit: g/cm^2
	Argument: 
		-dtogas: dust to gas ratio.
		-fraction: fraction of grain mass.
		-gas_surf_dens0: surface density of gas at reference radius.
	"""	
    def single_surf_dens0(self, dtogas, fraction, gas_surf_dens0):
	    return dtogas*fraction*gas_surf_dens0
	
	
    """ M)
	Equation (49). averaged grain size of the ith interval. Unit: cm
	Argument:
		-d_exp: size distribution exponent.
	"""
    def average(self, d_exp): 
	    return (((1-d_exp)/(4-d_exp))*((self.maxf**(4-d_exp) - self.minf**(4-d_exp))/(\
	    self.maxf**(1-d_exp) - self.minf**(1-d_exp))))**(1./3.) 
		

    """ N)
	Not written yet. mass of single grain.
	Argument:
		-d_exp: size distribution exponent.
	"""
    def single_mass(self, rho_m, size):
        return ((4.*np.pi)/3.)*rho_m*size**3
    
	##METHODS UNDER CONSTRUCTION
    """ O)
	q_ext
	"""
    def q_ext(self, A, lambd, size_av, q_c):
	    lambda_c = 2*np.pi*size_av
	    #print lambda_c
	    if lambd <= np.pi*size_av:
		    return A  #regime A
	    elif np.pi*size_av < lambd < 2*np.pi*size_av:
		    return q_c*(lambd/lambda_c)**(np.log10(q_c)/np.log10(2)) #regime B
	    elif lambd >= 2*np.pi*size_av:
		    return q_c*(lambd/lambda_c)**(-2) #regime C
		
		
    """ P)
	Av: pour une altitude pour l'instant.
	"""
    def Av(self, n_d, Q_ext, size_av, z):
        return 1.086*np.pi*n_d*Q_ext*z*size_av**2


#_________________________________________________#
#                      MAIN                       #
#_________________________________________________#
if __name__ == "__main__":
    #--initalization--
    i = j = total_dust_mass = total_dust_mass_out = 0
    total_dust_surdens0 = dust_surdens0 = fraction = dust_surdens = 0
    rho_d = rho_d_sum = 0
    ts_t = a_t = tot = 0
    lambd = 5e-5 #in "cm"

    nb_sizes = int(source_value(nb_sizes)) #8
    nb_points = int(source_value(nb_points)) #64
    print("number of sizes = %d" % nb_sizes)

    #--call Gas object--
    gas = Gas(source_value(Tmidplan_ref), source_value(ref_radius), source_value(star_mass)) 

    #--get C--
    grain_min = Grains(0, source_value(a_min), source_value(a_max), nb_sizes)
    grain_max = Grains(nb_sizes-1, source_value(a_min), source_value(a_max), nb_sizes)
    amin = grain_min.average(source_value(d_exp))
    amax = grain_max.average(source_value(d_exp))
    C_in = grain_min.C(source_value(dtogas_in), source_value(rho_m), source_value(d_exp), amax, amin)
    C_out = grain_min.C(source_value(dtogas_out), source_value(rho_m), source_value(d_exp), amax, amin)
    print("C_in = %.3E and C_out = %.3E" % (C_in, C_out))
    
	
    #--loop to get the extinction efficiency-- 
    for j in range(0, nb_sizes, 1):
	    grain = Grains(j, source_value(a_min), source_value(a_max), nb_sizes)
	    size_av = grain.average(source_value(d_exp)) 
	    grain.q_ext(1, lambd, size_av, 4)

    #--get the transition size in microns at the cut radius--
    grain0 = Grains(0, source_value(a_min), source_value(a_max), nb_sizes)
    ts_t = grain0.Ts_t(source_value(alpha), source_value(schmidt_number), source_value(s_x))
    a_t = grain0.a_t(ts_t, source_value(surface_density_ref), source_value(rho_m), source_value(cut_radius), source_value(ref_radius), source_value(p_exp))
    print("transition grain size: %E microns" % (a_t*1e4))
	
    #--get the gas scale heights in AU at reference radius--
    Href = gas.H_ref
    Hatmos_ref = gas.Hatmos_ref(source_value(Tatmos_ref), source_value(ref_radius), source_value(star_mass))
    print("gas scale height at reference radius: %.2f au" % (Href/cst.autocm))
	
    #--loop to retrieve the total dust mass-- 
    for i in range(0, nb_sizes, 1):
        grain = Grains(i, source_value(a_min), source_value(a_max), nb_sizes)
        total_dust_mass += grain.mass_density(source_value(d_exp), source_value(rho_m), C_in) # the total mass is 100% of the mass so this is used to get the mass fraction of each grain species.
    print(total_dust_mass)
    for size in range(0, nb_sizes, 1):
        grain = Grains(size, source_value(a_min), source_value(a_max), nb_sizes) # make sure to write this before anything in the loop.
        size_av = grain.average(source_value(d_exp))
        if size_av <= a_t:
            total_dust_mass_out += grain.mass_density(source_value(d_exp), source_value(rho_m), C_out) # the total mass is 100% of the mass so this is used to get the mass fraction of each grain species for outer radii.

    print(total_dust_mass_out)

    #--loop to retrieve the total grain surface density at reference radius--
    for size in range(0, nb_sizes, 1):
        grain = Grains(size, source_value(a_min), source_value(a_max), nb_sizes)
        size_av = grain.average(source_value(d_exp))
        fraction = (grain.mass_density(source_value(d_exp), source_value(rho_m), C_in)/total_dust_mass) # multiply by 100.0 to get the fraction in percentage of total dust mass.
        total_dust_surdens0 += grain.single_surf_dens0(source_value(dtogas_in), fraction, source_value(surface_density_ref))
        tot += fraction*100 #total in percentage. The value MUST be 100.
        print("size: %.5f microns -- fraction: %.5f -- cumul: %.5f" % (size_av*1e+4, fraction*100.0, tot))
    print("total fraction: %f" % tot)
	
    if int(source_value(settling)) == 0:
        print("settling: no")
    elif int(source_value(settling)) == 1:
        print("settling: yes")

   
    with open("radii.in","r") as radii:
        for line in radii:
            if line.strip():
                # re-initialize the quantity for each radius.
                n_H = n_d = H_d = rho_d = rho_g = A_V = nb = nb2 = 0 
                av = av2 = z_b = []
                dz = [0]

                
                parts = line.split()
                radius = parts[0]
                sigma_g = gas.surdens(source_value(surface_density_ref), float(radius), source_value(ref_radius), source_value(p_exp))
                sigma_g0 = gas.surdens(source_value(surface_density_ref), source_value(ref_radius), source_value(ref_radius), source_value(p_exp))
                H = gas.H(source_value(q_exp), float(radius), source_value(ref_radius)) #. divide by cst.autocm if "AU" instead of "cm".
                temp_mid = gas.temperature(source_value(Tmidplan_ref), float(radius), source_value(ref_radius), source_value(q_exp))
                temp_atmos = gas.temp_atmos(source_value(Tatmos_ref), float(radius), source_value(ref_radius), source_value(q_exp))
                mid = gas.dens_mid(sigma_g, H) 
                omega2 = gas.omega2(float(radius), source_value(star_mass))
                UV = gas.uv_factor(source_value(uv_ref), source_value(ref_radius), float(radius), H)
                print("inner radius = %d AU -- Hgas = %.5E AU -- surface density = %E -- Tmidplane = %.5E K -- Tatmos = %.5E K" % (int(radius), H/cst.autocm, sigma_g, temp_mid, temp_atmos))
 

                #----------------------------------------
                #--    creating the radius folders     --
                #----------------------------------------
                
                #--check if the folder exists, if not creates it--
                if not os.path.exists("MODEL/%sAU" % radius):
                    os.makedirs("MODEL/%sAU" % radius)


                #----------------------------------------
                #-- creating the "1D_static.dat" files --
                #----------------------------------------

                ##################### Get z, nH, Tg, Av, diff_coef, Td, 1/abund, AV_NH, size ##################

                z = gas.altitude2(nb_points, H) #!! in cm
                vertical_temperature = gas.vertical_temperature(temp_mid, temp_atmos, z, z[0], source_value(sigma_T))
                nH = gas.n_H2(mid, vertical_temperature, nb_points, omega2, z)
                dtog_up = gas.he_ab(source_value(dtogas_up), source_value(initial_He ))
                dtogas_mid = gas.dtogas_mid(nb_points, source_value(transition_altitude), source_value(dtogas_in), dtog_up, z, nH)
                grain_r, gtodm, AV_NH, Av_z = gas.gtodm(nb_points, H, nH, source_value(transition_altitude), source_value(rho_m), source_value(big_grains), source_value(small_grains), dtog_up, source_value(nH_to_AV_conversion), dtogas_mid, z)
                diff_coef = Td = np.zeros(64)
                

                ####################### Get the Av from dust ######################

                dz = gas.diff(z, nb_points)
                #___________________R <= Rc_________________
                if float(radius) <= source_value(cut_radius):
                    dtogas = source_value(dtogas_in)
                    for spatial_pt in range(0, nb_points,1):
                        a_v = 0
                        #write the Av column
                        for j in range(0, nb_sizes, 1):
                            grain = Grains(j, source_value(a_min), source_value(a_max), nb_sizes)
                            size_av = grain.average(source_value(d_exp))
                            fraction = grain.mass_density(source_value(d_exp), source_value(rho_m), C_in)/total_dust_mass
                            dust_surdens0 = grain.single_surf_dens0(dtogas, fraction, sigma_g0)
                            dust_surdens = grain.single_surf_dens(dust_surdens0, float(radius), source_value(ref_radius), source_value(p_exp))
                            single_mass = grain.single_mass(source_value(rho_m), size_av)
                            Ts0 = grain.stopping_time0(size_av, source_value(rho_m), sigma_g)
                            H_d = grain.H_d(H, Ts0, source_value(schmidt_number), source_value(alpha))
                            n_d = grain.rho_d(dust_surdens, H_d, z)/single_mass
                            q_ext = grain.q_ext(1, lambd, size_av, 4)
                            a_v = a_v + grain.Av(n_d[spatial_pt], q_ext, size_av, dz[spatial_pt])
                        A_V = A_V + a_v
                        av.append(A_V)
                        #print A_V
                    Av = np.asarray(av)
                    #print av2

                #___________________R > Rc__________________	
                elif float(radius) > source_value(cut_radius):
                    dtogas = source_value(dtogas_in) #sets the outer dust to gas ratio.		
                    for spatial_pt in range(0, nb_points,1):
                        a_v = 0
                        #write the Av column
                        for j in range(0, nb_sizes, 1):
                            grain = Grains(j, source_value(a_min), source_value(a_max), nb_sizes)
                            size_av = grain.average(source_value(d_exp))
                            if size_av <= a_t:
                                fraction = grain.mass_density(source_value(d_exp), source_value(rho_m), C_in)/total_dust_mass
                                dust_surdens0 = grain.single_surf_dens0(dtogas, fraction, source_value(surface_density_ref))
                                dust_surdens = grain.single_surf_dens(dust_surdens0, float(radius), source_value(ref_radius), source_value(p_exp))
                                single_mass = grain.single_mass(source_value(rho_m), size_av)
                                Ts0 = grain.stopping_time0(size_av, source_value(rho_m), sigma_g)
                                H_d = grain.H_d(H, Ts0, source_value(schmidt_number), source_value(alpha))
                                n_d = grain.rho_d(dust_surdens, H_d, z)/single_mass
                                q_ext = grain.q_ext(1, lambd, size_av, 4)
                                a_v = a_v + grain.Av(n_d[spatial_pt], q_ext, size_av, dz[spatial_pt])
                        A_V = A_V + a_v
                        av.append(A_V)
                    Av = np.asarray(av)
                    #print Av
                
                static_array = np.stack((z/cst.autocm, nH, vertical_temperature, Av_z, diff_coef, Td, gtodm, AV_NH, grain_r ), axis=-1)
                header_static = 'Reference radius AU = %.3E\n\
Midplan temperature (K) at the reference radius AU = %.3E\n\
Atmosphere (4H) temperature (K) at the reference radius AU = %.3E\n\
Radius where the structure is computed (AU) = %.3E\n\
Midplan temperature (K) at the requested radius = %.3E\n\
Atmosphere temperature (K) at the requested radius = %.3E\n\
Mass of the central object (g) = %.3E\n\
UV field coming from the star at the reference radius AU (unit: mean ISM field) = %.3E\n\
UV field coming from the star at the requested radius = %.3E\n\
Vertical resolution =  %d\n\
Scale Height H in au = %.3E\n\
Distance [AU] ; H Gas density [part/cm^3] ; Gas Temperature [K] ; Visual Extinction [mag] ; Diffusion coefficient [cm^2/s]; Dust Temperature [K]; 1/abundance of grains ; AV/NH conversion factor ; radius of grains (cm)'\
    % (source_value(ref_radius), source_value(Tmidplan_ref), source_value(Tatmos_ref), float(radius), vertical_temperature[63], vertical_temperature[0], cst.Msun*source_value(star_mass), source_value(uv_ref), UV, nb_points, H/cst.autocm)
                np.savetxt('MODEL/%sAU/1D_static.dat' % radius, static_array, fmt='%.5E', delimiter='   ', newline='\n', header=header_static , comments='! ', encoding=None)
                

                #_________________Remove blank lines, very important!________________
                lines = [line.strip() for line in fileinput.FileInput('MODEL/%sAU/1D_static.dat' % radius) ]

                #----------------------------------------
                #-- updating the "parameters.in" files --
                #----------------------------------------

                #_________________R <= Rc________________
                if float(radius) <= source_value(cut_radius):
                    sed(NBR_GRAINS, str(nb_sizes), "parameters.in", "MODEL/%sAU/parameters.in" % radius)
                    sed(UV_FLUX, "%.3E" %UV, "MODEL/%sAU/parameters.in" % radius)
                

                #__________________R > Rc________________
                elif float(radius) > source_value(cut_radius):
                    for size in range(0, nb_sizes, 1):
                        grain = Grains(size, source_value(a_min), source_value(a_max), nb_sizes) # make sure to write this before anything in the loop.
                        size_av = grain.average(source_value(d_exp))
                        if size_av <= a_t:
                            nb+=1	
                    sed(NBR_GRAINS, str(nb), "parameters.in", "MODEL/%sAU/parameters.in" % radius)
                    sed(UV_FLUX, "%.3E" %UV, "MODEL/%sAU/parameters.in" % radius)

                
                #----------------------------------------
                #-- creating the "1D_grain_size" files --
                #----------------------------------------

                with open("MODEL/%sAU/1D_grain_sizes.in" % radius, "w+") as abundance:
                    header_grain = "! grain-radius [cm]                     \
                    1/abundance of grains           \
                    grain temp     \
                    CR-peak-Temperaturegrain[K]      \
                    spatial point\
                    "
                    abundance.write(header_grain)
                    abundance.write("\n\n")


                    #__________R <= Rc_________
                    if float(radius) <= source_value(cut_radius):
                        dtogas = source_value(dtogas_in) #sets the inner dust to gas ratio.
                        for spatial_pt in range(0, nb_points,1):
                            
                            #write the grain sizes 
                            for size in range(0, nb_sizes, 1):
                                grain = Grains(size, source_value(a_min), source_value(a_max), nb_sizes) # make sure to write this before anything in the loop.
                                size_av = grain.average(source_value(d_exp))
                                abundance.write("%.4E " % size_av)
                            abundance.write("    ")
                            
                            
                            if int(source_value(settling)) == 0.000e+00:
                                #write the grain densities
                                for size in range(0, nb_sizes, 1):	
                                    grain = Grains(size, source_value(a_min), source_value(a_max), nb_sizes)
                                    size_av = grain.average(source_value(d_exp))
                                    fraction = grain.mass_density(source_value(d_exp), source_value(rho_m), C_in)/total_dust_mass
                                    single_mass = grain.single_mass(source_value(rho_m), size_av)
                                    n_d = np.asarray((fraction*dtogas*nH*cst.mu*cst.amu)/single_mass)
                                    #abundance = n_d/n_H
                                    #print "radius = %s, H_d = %e, sigma_d = %e" % (radius, H_d, dust_surdens)
                                    #if n_d[spatial_pt] <= 1.000000E-70:
                                        #abundance.write("1.000000E+80 ")
                                    #else:
                                    n_d[spatial_pt] = max(n_d[spatial_pt], 1e-45)
                                    abundance.write("%E " % np.true_divide(nH[spatial_pt], n_d[spatial_pt]))		
                                abundance.write("\t")
                            
                            elif int(source_value(settling)) == 1.000e+00:
                                #write the grain densities
                                for size in range(0, nb_sizes, 1):	
                                    grain = Grains(size, source_value(a_min), source_value(a_max), nb_sizes)
                                    size_av = grain.average(source_value(d_exp))
                                    fraction = grain.mass_density(source_value(d_exp), source_value(rho_m), C_in)/total_dust_mass
                                    dust_surdens0 = grain.single_surf_dens0(dtogas, fraction, source_value(surface_density_ref))
                                    dust_surdens = grain.single_surf_dens(dust_surdens0, float(radius), source_value(ref_radius), source_value(p_exp))
                                    single_mass = grain.single_mass(source_value(rho_m), size_av)
                                    Ts0 = grain.stopping_time0(size_av, source_value(rho_m), sigma_g)
                                    H_d = grain.H_d(H, Ts0, source_value(schmidt_number), source_value(alpha))
                                    n_d = np.asarray(grain.rho_d(dust_surdens, H_d, z)/single_mass)
                                    #abundance = n_d/n_H
                                    #print "radius = %s, H_d = %e, sigma_d = %e" % (radius, H_d, dust_surdens)
                                    #if n_d[spatial_pt] <= 1.000000E-70:
                                        #abundance.write("1.000000E+80 ")
                                    #else:
                                    n_d[spatial_pt] = max(n_d[spatial_pt], 1e-45)
                                    abundance.write("%E " % np.true_divide(nH[spatial_pt], n_d[spatial_pt]))		
                                abundance.write("\t")
                            
                            #write the grain temperature
                            temperature = np.loadtxt('../SOURCE/temperature/T%s.txt' % radius, unpack=True)
                            for size in range(0, nb_sizes, 1):
                                abundance.write("%e " % temperature[size,spatial_pt]) #column, row
                            abundance.write("\t")
                            
                            #write the grain temperature 
                            #for size in range(0, nb_sizes, 1):
                                #abundance.write("10 ")
                            #abundance.write("\t")
                            
                            #write CR-peak-temperature grain 
                            for size in range(0, nb_sizes, 1):
                                abundance.write("73 ")
                            abundance.write("            !         %s" % (spatial_pt+1))
                            abundance.write("\n")
                
                    #__________R > Rc_________	
                    elif float(radius) > source_value(cut_radius):
                        dtogas = source_value(dtogas_in) #sets the outer dust to gas ratio.
                        for spatial_pt in range(0, nb_points,1):
                            
                            #write the grain sizes
                            for size in range(0, nb_sizes, 1):
                                grain = Grains(size, source_value(a_min), source_value(a_max), nb_sizes) # make sure to write this before anything in the loop.
                                size_av = grain.average(source_value(d_exp))
                                if size_av <= a_t:
                                    abundance.write("%.2E " % size_av)	
                            abundance.write("     ")
                            
                            if int(source_value(settling)) == 0.000e+00:
                                #write the grain densities
                                for size in range(0, nb_sizes, 1):
                                    grain = Grains(size, source_value(a_min), source_value(a_max), nb_sizes) # make sure to write this before anything in the loop.
                                    size_av = grain.average(source_value(d_exp))
                                
                                    if size_av <= a_t:
                                        fraction = grain.mass_density(source_value(d_exp), source_value(rho_m), C_in)/total_dust_mass
                                        single_mass = grain.single_mass(source_value(rho_m), size_av)
                                        n_d = np.asarray((fraction*dtogas*nH*cst.mu*cst.amu)/single_mass)
                                        #abundance = n_d/n_H
                                        #if n_d[spatial_pt] <= 1.000000E-70:
                                            #abundance.write("1.000000E+80 ")
                                        #else:
                                        n_d[spatial_pt] = max(n_d[spatial_pt], 1e-45)
                                        abundance.write("%E " % np.true_divide(nH[spatial_pt], n_d[spatial_pt]))	
                                abundance.write("\t\t")
                                
                            elif int(source_value(settling)) == 1.000e+00:
                                #write the grain densities
                                for size in range(0, nb_sizes, 1):
                                    grain = Grains(size, source_value(a_min), source_value(a_max), nb_sizes) # make sure to write this before anything in the loop.
                                    size_av = grain.average(source_value(d_exp))
                                
                                    if size_av <= a_t:
                                        fraction = grain.mass_density(source_value(d_exp), source_value(rho_m), C_in)/total_dust_mass
                                        dust_surdens0 = grain.single_surf_dens0(dtogas, fraction, source_value(surface_density_ref))
                                        dust_surdens = grain.single_surf_dens(dust_surdens0, float(radius), source_value(ref_radius), source_value(p_exp))
                                        single_mass = grain.single_mass(source_value(rho_m), size_av)
                                        Ts0 = grain.stopping_time0(size_av, source_value(rho_m), sigma_g)
                                        H_d = grain.H_d(H, Ts0, source_value(schmidt_number), source_value(alpha))
                                        n_d = np.asarray(grain.rho_d(dust_surdens, H_d, z)/single_mass)
                                        #abundance = n_d/n_H
                                        #if n_d[spatial_pt] <= 1.000000E-70:
                                            #abundance.write("1.000000E+80 ")
                                        #else:
                                        n_d[spatial_pt] = max(n_d[spatial_pt], 1e-45)
                                        abundance.write("%E " % np.true_divide(nH[spatial_pt], n_d[spatial_pt]))	
                                abundance.write("\t\t")
                            
                            #write the grain temperature 
                            temperature = np.loadtxt('../SOURCE/temperature/T%s.txt' % radius, unpack=True)
                            for size in range(0, nb_sizes, 1):
                                grain = Grains(size, source_value(a_min), source_value(a_max), nb_sizes) # make sure to write this before anything in the loop.
                                size_av = grain.average(source_value(d_exp))
                                if size_av <= a_t:
                                    abundance.write("%e " % temperature[size,spatial_pt]) #column, row
                            abundance.write("\t")
                            
                            #write the grain temperature 
                            #for size in range(0, nb_sizes, 1):
                                #grain = Grains(size) # make sure to write this before anything in the loop.
                                #size_av = grain.average(source_value(d_exp))
                                #if size_av <= a_t:
                                    #abundance.write("10 ")
                            #abundance.write("\t")
                            
                            #write CR-peak-temperature grain 
                            for size in range(0, nb_sizes, 1):
                                grain = Grains(size, source_value(a_min), source_value(a_max), nb_sizes) # make sure to write this before anything in the loop.
                                size_av = grain.average(source_value(d_exp))
                                if size_av <= a_t:
                                    abundance.write("73 ")
                            abundance.write("            !         %s" % (spatial_pt+1))
                            abundance.write("\n")
                abundance.close()


                #---------------------------------------------
                #-- creating "dust_density" and sizes files --
                #---------------------------------------------
 
                #__________R <= Rc_________
                if float(radius) <= source_value(cut_radius):
                    dtogas = source_value(dtogas_in) #sets the inner dust to gas ratio.
                    n_d_array = np.array([])
                    name_list = np.array([])
                    size_list = np.array([])
                    #write the grain sizes 
                    for size in range(0, nb_sizes, 1):
                        grain = Grains(size, source_value(a_min), source_value(a_max), nb_sizes) # make sure to write this before anything in the loop.
                        size_av = grain.average(source_value(d_exp))
                        fraction = grain.mass_density(source_value(d_exp), source_value(rho_m), C_in)/total_dust_mass
                        single_mass = grain.single_mass(source_value(rho_m), size_av)
                        dust_surdens0 = grain.single_surf_dens0(dtogas, fraction, source_value(surface_density_ref))
                        dust_surdens = grain.single_surf_dens(dust_surdens0, float(radius), source_value(ref_radius), source_value(p_exp))
                        Ts0 = grain.stopping_time0(size_av, source_value(rho_m), sigma_g)
                        H_d = grain.H_d(H, Ts0, source_value(schmidt_number), source_value(alpha))
                        name_list = np.append(name_list, "size%d" % (size+1))
                        size_list = np.append(size_list, size_av)
                        if int(source_value(settling)) == 0.000e+00:
                            n_d = np.asarray((fraction*dtogas*nH*cst.mu*cst.amu)/single_mass)
                        elif int(source_value(settling)) == 1.000e+00:
                            n_d = np.asarray(grain.rho_d(dust_surdens, H_d, z)/single_mass)
                        n_d = np.maximum(n_d, 1e-45)
                        n_d_array = np.hstack((n_d_array, n_d))
                    size_list = np.vstack((name_list, size_list))
                    name_str = " ".join(name_list )	
                    n_d_array = np.reshape(n_d_array, (nb_sizes, nb_points))
                    np.savetxt('MODEL/%sAU/dust_density.dat' % radius, n_d_array.T, fmt='%.5E', delimiter=' ', newline='\n', header=name_str, comments='', encoding=None)
                    np.savetxt('MODEL/%sAU/sizes.dat' % radius, size_list.T, fmt="%s", delimiter=' ', newline='\n', encoding=None)

                #__________R > Rc_________	
                elif float(radius) > source_value(cut_radius):
                    dtogas = source_value(dtogas_out) #sets the outer dust to gas ratio.
                    n_d_array = np.array([])
                    name_list = np.array([])
                    size_list = np.array([])
                    #write the grain sizes
                    for size in range(0, nb_sizes, 1):
                        grain = Grains(size, source_value(a_min), source_value(a_max), nb_sizes) # make sure to write this before anything in the loop.
                        size_av = grain.average(source_value(d_exp))
                        if size_av <= a_t:
                            nb2+=1
                            fraction = grain.mass_density(source_value(d_exp), source_value(rho_m), C_in)/total_dust_mass_out
                            single_mass = grain.single_mass(source_value(rho_m), size_av)
                            dust_surdens0 = grain.single_surf_dens0(dtogas, fraction, source_value(surface_density_ref))
                            dust_surdens = grain.single_surf_dens(dust_surdens0, float(radius), source_value(ref_radius), source_value(p_exp))
                            Ts0 = grain.stopping_time0(size_av, source_value(rho_m), sigma_g)
                            H_d = grain.H_d(H, Ts0, source_value(schmidt_number), source_value(alpha))
                            name_list = np.append(name_list, "size%d" % (size+1))
                            size_list = np.append(size_list, size_av)
                            if int(source_value(settling)) == 0.000e+00:
                                n_d = np.asarray((fraction*dtogas*nH*cst.mu*cst.amu)/single_mass)
                            elif int(source_value(settling)) == 1.000e+00:
                                n_d = np.asarray(grain.rho_d(dust_surdens, H_d, z)/single_mass)
                            n_d = np.maximum(n_d, 1e-45)
                            n_d_array = np.hstack((n_d_array, n_d))
                    size_list = np.vstack((name_list, size_list))
                    name_str = " ".join(name_list )	
                    n_d_array = np.reshape(n_d_array, (nb2, nb_points))
                    np.savetxt('MODEL/%sAU/dust_density.dat' % radius, n_d_array.T, fmt='%.5E', delimiter=' ', newline='\n', header=name_str, comments='', encoding=None)
                    np.savetxt('MODEL/%sAU/sizes.dat' % radius, size_list.T, fmt="%s", delimiter=' ', newline='\n', encoding=None)

                #----------------------------------------
                #--     transfer the initial files     --
                #----------------------------------------

                for filename in glob.glob(os.path.join("../SOURCE/initial_files/", "*.*")):
                    copy(filename, "MODEL/%sAU/" % radius)


    """
    static = np.loadtxt("100AU.txt", comments='!', skiprows=0)
    av_naut = static[:,3]
    print av_naut

    plt.xlabel('Av', fontsize = 17)
    plt.ylabel('elevation z (AU)', fontsize = 17)
    plt.ylim(0,4)
    plt.semilogy(z/(H), nH, linewidth=1, color="black", label="T from Williams and Best (2014)")
    plt.semilogx(av2, z/(H), linewidth=1, color="black", linestyle='-', label="T from Williams and Best (2014)")
    plt.semilogx(av_naut, z/(H), linewidth=1, color="black", linestyle='--', label="T from Williams and Best (2014)")
    plt.legend()
    plt.show()
    """