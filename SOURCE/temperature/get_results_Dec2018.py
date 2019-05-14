# Code for extracting temperature and density distribution from model.fits-file


import argparse
import numpy as np
from astropy.io import fits

#=======================================================================
#required model and grid parameters, DO NOT CHANGE!
n_species = 16
R_IN = 1				#inner disk radius
R_OU = 300				#outer disk radius
N_A = 300				#number of cells in r direction
sf_A = 1.03				#stepwidthfaktor r direction
N_B_ou = 55
N_B = 81 +2*N_B_ou 

#=======================================================================

def in_model_space(r, z):
	"""
	Test if the given position is inside model space
	"""
	
	radius = np.sqrt(r**2 + z**2)
	ims = (radius >= R_IN) & (radius <= R_OU)
	
	flag = False
	for num, is_in_model_space in enumerate(ims):
		if not is_in_model_space:
			print('coordinate (' + str(r[num]) + '/' + str(z[num]) + ') not in model space')
			flag = True
	if flag:
		print('ERROR: some coordinates are outside the model space (sperical model space with Rin=' 
		+ str(R_IN) + ' au, Rou=' + str(R_OU) + ' au)')
		exit()
    
def calc_cell_number(r, z):		
	"""
	r: array of radius values, z:array of corresponding height values
	"""
	#Check if all coordinates are in model space:
	in_model_space(r, z)	
	#calc cell boundaries in r direction:
	a_bounds = np.zeros(N_A+1, dtype=np.float64)
	a1=(R_OU-R_IN)*(sf_A-1)/(sf_A**N_A-1)
	for i in range(N_A+1):
		a_bounds[i]=R_IN+a1*(sf_A**i-1)/(sf_A-1)
	a_midpoints = (a_bounds[1:] + a_bounds[:-1]) / 2
		
	#calc cell boundaries in z direction (theta coordinate (-PI/2,...,+PI/2) linear in [-0.3 rad, 0.3 rad] log else):
	b_bounds_mi = np.linspace(-0.3, 0.3, N_B + 1 - 2*N_B_ou, dtype=np.float64)
	b_bounds_ou = np.logspace(np.log10(0.3), np.log10(np.pi/2), N_B_ou+1, dtype=np.float64)
	
	b_bounds=np.zeros(N_B + 1, dtype=np.float64)
	b_bounds[0:N_B_ou]=-np.flipud(b_bounds_ou[1:N_B_ou+1])
	b_bounds[N_B_ou:N_B + 1-N_B_ou]=b_bounds_mi
	b_bounds[N_B -(N_B_ou-1): N_B+1]=b_bounds_ou[1:N_B_ou+1]
	b_midpoints = (b_bounds[1:] + b_bounds[:-1]) / 2
	#find cell number:	
	n_cell = []
	midpoint_r = []
	midpoint_z =[]
	radius = np.sqrt(r**2 + z**2)
	theta = np.arctan2(z,r)
	for rval, theta_val in zip(radius, theta):
		if rval>=a_bounds[-1]:
			n_r = N_A-1
		else:
			n_r = np.amin(np.argwhere(np.array(a_bounds)>rval))
			
		if theta_val>=b_bounds[-1]:
			n_theta = N_B-1
		else:
			n_theta = np.amin(np.argwhere(np.array(b_bounds)>theta_val))
		midpoint_r.append(a_midpoints[n_r-1] * np.cos(b_midpoints[n_theta-1]))
		midpoint_z.append(a_midpoints[n_r-1] * np.sin(b_midpoints[n_theta-1]))
		n_cell.append((n_r-1)*N_B+n_theta)
	n_cell = np.array(n_cell)
	midpoint_r = np.array(midpoint_r)
	midpoint_z = np.array(midpoint_z)
	return n_cell, midpoint_r, midpoint_z
	
def read_temp_den_from_modelfits(fitsfile, n_cell, r, z, midpoint_r, midpoint_z):
	"""
	read density and temperature values for given coordinates from model file
	"""
	hdu_list = fits.open(fitsfile)
	disk_data = hdu_list[0].data
	hdu_list.close()
	
	density_data = disk_data[n_cell-1, 3:3+n_species]
	temp_data = disk_data[n_cell-1, 3+n_species+4:3+n_species+4+n_species] 
	if len(n_cell)==1:
		for i in range(n_species):
			print('grain size interval ' + str(i+1) + ' : den = ' + str(density_data[0,i]) + ' m^-3, T = ' + str(temp_data[0,i]) + ' K')
		print("cell midpoint: r = " + str(midpoint_r[0]) + " au, z = " + str(midpoint_z[0]) + " au")
	else:
		den_array = np.zeros((len(n_cell), n_species + 4))
		den_array[:,0]=r
		den_array[:,1]=z
		den_array[:,2]=midpoint_r
		den_array[:,3]=midpoint_z
		den_array[:,4:] = density_data
		np.savetxt('density.txt', den_array, header="r[au]    z[au]  nearest cell r[au] nearest cell z[au]  number density of grain size interval 1 [m^-3]    ...    number density of grain size interval 16 [m^-3]")
		temp_array = np.zeros((len(n_cell), n_species + 4))
		temp_array[:,0]=r
		temp_array[:,1]=z
		temp_array[:,2]=midpoint_r
		temp_array[:,3]=midpoint_z
		temp_array[:,4:]=temp_data
		np.savetxt('temperature.txt', temp_array, header="r[au]    z[au]	nearest cell r[au] nearest cell z[au]    temperature of grain size interval 1 [K]    ...    temperature of grain size interval 16 [K]")
		
def read_coordinates_from_file(coordinates):
	"""
	read requested coordinates from file
	"""
	coords = np.loadtxt(coordinates)
	r = coords[:,0]
	z = coords[:,1]
	
	return r, z


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Returns dust number density [m^-3] and temperature [K] values for given coordinates.')
	
	parser.add_argument('-f', '--fitsfile', type=str,
						default='model.fits',
						help='Enter the fits-file name, from which density and temperature values are to be extracted. ')
						
	parser.add_argument('-c', '--coordfile', dest='coordinates', type=str,
						help='Enter file with cartesian coordinates, 1 column: radius r [au], 2. column height z [au]')
						
	parser.add_argument('--coord', dest='coordinate', type=float, nargs=2,
						help='Enter cartesian coordinates radius r [au] and height z [au]')
	
	args = parser.parse_args()

	if args.coordinates is not None:
		print("Reading density and temperature values for coordinates from " + args.coordinates)
		r, z = read_coordinates_from_file(args.coordinates)
	elif args.coordinate is not None:
		print("Reading density and temperature values for (" + str(args.coordinate[0]) + 'au/' + str(args.coordinate[1]) + 'au)' )
		r = np.array([args.coordinate[0]])
		z = np.array([args.coordinate[1]])
	else:
		print("Please give coordinates, see help for details.")
		
	n_cell, midpoint_r, midpoint_z = calc_cell_number(r, z)
	read_temp_den_from_modelfits(args.fitsfile, n_cell, r, z, midpoint_r, midpoint_z)
	
