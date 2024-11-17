import matplotlib.pyplot as plt
import numpy as np
import os
import scipy as sp
import scipy.optimize as opt

from lmfit.models import GaussianModel

filelist = ["alpha_Fe.TKA", 
			"Fe2O3.TKA", 
			"K4Fe_CN_6.TKA", 
			"LiFePO4.TKA", 
			"stainless.TKA"]

missing_files = 0
for f in filelist:
	if not os.path.isfile(f):
		print("Missing file: ", f)
	
if missing_files: exit()

#calibrate = [90, 180, 55, 65, 135]
calibrate = [136, 89, 46, 66, 108]
exp_peaks = [12, 12, 2, 4, 2]

def gen_model(num_peaks, guess_centers, guess_amps):
	if num_peaks == 2:
		model = (GaussianModel(prefix='p0_')
				+GaussianModel(prefix='p1_'))

		params = model.make_params(p0_center=guess_centers[0], p0_amplitude=guess_amps[0], p0_sigma=0.005,
									p1_center=guess_centers[1], p1_amplitude=guess_amps[1], p1_sigma=0.005)
		return model, params
	elif num_peaks == 4:
		model = (GaussianModel(prefix='p0_')
				+GaussianModel(prefix='p1_')
				+GaussianModel(prefix='p2_')
				+GaussianModel(prefix='p3_'))

		params = model.make_params(p0_center=guess_centers[0], p0_amplitude=guess_amps[0], p0_sigma=0.005,
									p1_center=guess_centers[1], p1_amplitude=guess_amps[1], p1_sigma=0.005,
									p2_center=guess_centers[2], p2_amplitude=guess_amps[2], p2_sigma=0.005,
									p3_center=guess_centers[3], p3_amplitude=guess_amps[3], p3_sigma=0.005)
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
		params = model.make_params(p0_center=guess_centers[0], p0_amplitude=guess_amps[0], p0_sigma=0.005,
									p1_center=guess_centers[1], p1_amplitude=guess_amps[1], p1_sigma=0.005,
									p2_center=guess_centers[2], p2_amplitude=guess_amps[2], p2_sigma=0.005,
									p3_center=guess_centers[3], p3_amplitude=guess_amps[3], p3_sigma=0.005,
									p4_center=guess_centers[4], p4_amplitude=guess_amps[4], p4_sigma=0.005,
									p5_center=guess_centers[5], p5_amplitude=guess_amps[5], p5_sigma=0.005,
									p6_center=guess_centers[6], p6_amplitude=guess_amps[6], p6_sigma=0.005,
									p7_center=guess_centers[7], p7_amplitude=guess_amps[7], p7_sigma=0.005,
									p8_center=guess_centers[8], p8_amplitude=guess_amps[8], p8_sigma=0.005,
									p9_center=guess_centers[9], p9_amplitude=guess_amps[9], p9_sigma=0.005,
									pa_center=guess_centers[10], pa_amplitude=guess_amps[10], pa_sigma=0.005,
									pb_center=guess_centers[11], pb_amplitude=guess_amps[11], pb_sigma=0.005)	
		return model, params

def sin_fit(t, c, A, w, phi):
	return c+A*np.sin(w*t+phi)
	
def exact_fit(t, a, b, c, w, phi):
	return a - b + np.sqrt(b**2-(c*np.cos(w*t+phi))**2)-c*np.sin(w*t+phi)

def gauss_fit(x, a, mu, sigma):
	return a*np.exp(-(x-mu)*(x-mu)/(2*sigma*sigma))

for index, f in enumerate(filelist):
	plt.figure(index)
	counts = np.loadtxt(f)[2:14262]	# so that data is 14260 = 2*2*5*23*31, can rebin into nice bin sizes
	bins = np.arange(len(counts))*70e-6+35e-6
	
	plt.subplot(3, 1, 1)
	plt.stairs(counts, np.arange(len(counts)+1)*70e-6, color='k')
	plt.xlim(left=0,right=1)
	
	'''
	popt, pcov = opt.curve_fit(sin_fit, bins, counts, 
								p0=[np.mean(counts), (max(counts)-min(counts))/2, 2*np.pi, 0])
	
	smooth = sp.signal.savgol_filter(sin_fit(bins, *popt)-counts, 51, 3)
	'''
	smooth = sp.signal.savgol_filter(-counts, 51, 3)
	peaks, _ = sp.signal.find_peaks(smooth, distance=100, prominence=calibrate[index])
	
	counts_filtered = np.copy(counts)
	bins_filtered = np.copy(bins)
	w = 100
	print(len(peaks))

	for peak in peaks:
		counts_filtered[(peak-w):(peak+w)] = np.nan
		bins_filtered[(peak-w):(peak+w)] = np.nan
	
	plt.subplot(3, 1, 2)
	plt.stairs(counts_filtered, np.arange(len(counts)+1)*70e-6, color='k')
	
	counts_filtered = counts_filtered[~np.isnan(counts_filtered)]
	bins_filtered = bins_filtered[~np.isnan(bins_filtered)]
	
	popt, pcov = opt.curve_fit(sin_fit, bins_filtered, counts_filtered, 
								p0=[np.mean(counts_filtered), (max(counts_filtered)-min(counts_filtered))/2, 2*np.pi, 0])
	
	plt.plot(bins, sin_fit(bins, *popt), color='r')
	plt.xlim(left=0,right=1)
	
	plt.subplot(3, 1, 3)
	counts_cleaned = sin_fit(bins, *popt)-counts
	plt.stairs(counts_cleaned, np.arange(len(counts)+1)*70e-6)
	plt.xlim(left=0,right=1)
	
	model, params = gen_model(exp_peaks[index], bins[peaks], counts_cleaned[peaks])
	result = model.fit(counts_cleaned, params, x=bins)
	
	#print(result.fit_report())

	plt.plot(bins, result.best_fit)
	
plt.show()
