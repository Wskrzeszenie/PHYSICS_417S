import matplotlib.pyplot as plt
import numpy as np
import os
import scipy as sp
import scipy.optimize as opt

from lmfit.models import GaussianModel

filelist = ["alpha_Fe", 
			"Fe2O3", 
			"K4Fe_CN_6", 
			"LiFePO4_1", 
			"stainless",
			"FeC2O4_1",
			"7mu_Fe",
			"LiFePO4_2",
			"FeC2O4_2"]

exp_peaks = [12, 12, 2, 4, 2, 4, 12, 4, 4]
c = 299792458 # m/s

def file_check():
	missing_files = 0
	for f in filelist:
		if not os.path.isfile(f+".TKA"):
			print("Missing file: ", f)
		
	if missing_files: exit()
	
def energy_approx(v):
	E = 14.4e3
	return E*v/c
	
def v(t, w, phi, model=0):
	r = 2e-3
	x = w*t+phi
	if model == 0:
		return r*w*np.cos(x)
	else:
		return r*r*np.cos(x)*np.sin(x)/np.sqrt(c*c-r*r*np.cos(x)*np.cos(x))+r*np.cos(x)
	return 0

def sin_fit(t, c, A, w, phi):
	return c+A*np.sin(w*t+phi)
	
#def exact_fit(t, a, b, c, w, phi):
#	return a - b + np.sqrt(b**2-(c*np.cos(w*t+phi))**2)-c*np.sin(w*t+phi)

def gauss_fit(x, a, mu, sigma):
	return a*np.exp(-(x-mu)*(x-mu)/(2*sigma*sigma))
	
def calculate_energies(result, popt):
	result_params = list(reversed(result.best_values.keys()))
	#total_fit = np.zeros_like(bins)
	for i in range(len(result_params)//3):
		amp = result.best_values[result_params[int(i*3)]]
		mu = result.best_values[result_params[int(i*3+1)]]
		sigma = result.best_values[result_params[int(i*3+2)]]
		print(np.round(mu, 3), 1e8*energy_approx(v(mu, popt[2], popt[3])))
		#total_fit += gauss_fit(bins, amp, mu, sigma)
	
def gen_model(num_peaks, guess_centers, guess_amps):
	if num_peaks == 2:
		model = (GaussianModel(prefix='p0_')
				+GaussianModel(prefix='p1_'))

		params = model.make_params(p0_center=guess_centers[0], p0_amplitude=guess_amps[0], p0_sigma=0.001,
									p1_center=guess_centers[1], p1_amplitude=guess_amps[1], p1_sigma=0.001)
		return model, params
	elif num_peaks == 4:
		model = (GaussianModel(prefix='p0_')
				+GaussianModel(prefix='p1_')
				+GaussianModel(prefix='p2_')
				+GaussianModel(prefix='p3_'))

		params = model.make_params(p0_center=guess_centers[0], p0_amplitude=guess_amps[0], p0_sigma=0.001,
									p1_center=guess_centers[1], p1_amplitude=guess_amps[1], p1_sigma=0.001,
									p2_center=guess_centers[2], p2_amplitude=guess_amps[2], p2_sigma=0.001,
									p3_center=guess_centers[3], p3_amplitude=guess_amps[3], p3_sigma=0.001)
		return model, params
	elif num_peaks == 12:						
		model = (GaussianModel(prefix='p0_')
				+GaussianModel(prefix='p1_')
				+GaussianModel(prefix='p2_')
				+GaussianModel(prefix='p3_')
				+GaussianModel(prefix='p4_')
				+GaussianModel(prefix='p5_')
				+GaussianModel(prefix='p6_')
				+GaussianModel(prefix='p7_')
				+GaussianModel(prefix='p8_')
				+GaussianModel(prefix='p9_')
				+GaussianModel(prefix='pa_')
				+GaussianModel(prefix='pb_')) 
		params = model.make_params(p0_center=guess_centers[0], p0_amplitude=guess_amps[0], p0_sigma=0.001,
									p1_center=guess_centers[1], p1_amplitude=guess_amps[1], p1_sigma=0.001,
									p2_center=guess_centers[2], p2_amplitude=guess_amps[2], p2_sigma=0.001,
									p3_center=guess_centers[3], p3_amplitude=guess_amps[3], p3_sigma=0.001,
									p4_center=guess_centers[4], p4_amplitude=guess_amps[4], p4_sigma=0.001,
									p5_center=guess_centers[5], p5_amplitude=guess_amps[5], p5_sigma=0.001,
									p6_center=guess_centers[6], p6_amplitude=guess_amps[6], p6_sigma=0.001,
									p7_center=guess_centers[7], p7_amplitude=guess_amps[7], p7_sigma=0.001,
									p8_center=guess_centers[8], p8_amplitude=guess_amps[8], p8_sigma=0.001,
									p9_center=guess_centers[9], p9_amplitude=guess_amps[9], p9_sigma=0.001,
									pa_center=guess_centers[10], pa_amplitude=guess_amps[10], pa_sigma=0.001,
									pb_center=guess_centers[11], pb_amplitude=guess_amps[11], pb_sigma=0.001)	
		return model, params
