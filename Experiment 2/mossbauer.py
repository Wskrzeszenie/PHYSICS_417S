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
			"stainless.TKA",
			"FeC2O4.TKA",
			"7mu_Fe.TKA"]

missing_files = 0
for f in filelist:
	if not os.path.isfile(f):
		print("Missing file: ", f)
	
if missing_files: exit()

#calibrate = [90, 180, 55, 65, 135]
#calibrate = [35, 25, 35, 35, 35, 35, 20]
calibrate = [12, 16, 16, 14, 7, 20, 16]
exp_peaks = [12, 12, 2, 4, 2, 4, 12]

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

def sin_fit(t, c, A, w, phi):
	return c+A*np.sin(w*t+phi)
	
def exact_fit(t, a, b, c, w, phi):
	return a - b + np.sqrt(b**2-(c*np.cos(w*t+phi))**2)-c*np.sin(w*t+phi)

def gauss_fit(x, a, mu, sigma):
	return a*np.exp(-(x-mu)*(x-mu)/(2*sigma*sigma))/(sigma*np.sqrt(2*np.pi))

plt.figure()

for index, f in enumerate(filelist):
	#plt.subplot(len(filelist),1,index+1)
	plt.suptitle(f)
	plt.figure(index)
	counts = np.loadtxt(f)[2:14262]	# so that data is 14260 = 2*2*5*23*31, can rebin into nice bin sizes
	bins = np.arange(len(counts))*70e-6+35e-6
	
	plt.subplot(3, 1, 1)
	plt.stairs(counts, np.arange(len(counts)+1)*70e-6, color='k')
	plt.xlim(left=0,right=1)
	
	n = 460
	bins_r = np.arange(n)*14620/n*70e-6+14620/n*35e-6
	counts_r = np.zeros(n)
	for i in range(14260):
		counts_r[i//(14260//n)] += counts[i]
	
	popt, pcov = opt.curve_fit(sin_fit, bins_r, counts_r, 
								p0=[np.mean(counts_r), (max(counts_r)-min(counts_r))/2, 2*np.pi, 0])
	
	counts_guess = sin_fit(bins_r, *popt)-counts_r
	peaks, _ = sp.signal.find_peaks(counts_guess/max(counts_guess)*100, distance = 5, prominence=calibrate[index])
	
	counts_rf = np.copy(counts_r)
	bins_rf = np.copy(bins_r)
	w = 3
	print(len(peaks))
	for peak in peaks:
		counts_rf[(peak-w):(peak+w)] = np.nan
		bins_rf[(peak-w):(peak+w)] = np.nan
		
	plt.subplot(3, 1, 2)
	plt.stairs(counts_rf, np.arange(n+1)*14620/n*70e-6, color='k')
	
	counts_rf = counts_rf[~np.isnan(counts_rf)]
	bins_rf = bins_rf[~np.isnan(bins_rf)]
	
	popt, pcov = opt.curve_fit(sin_fit, bins_rf, counts_rf, 
								p0=[np.mean(counts_rf), (max(counts_rf)-min(counts_rf))/2, 2*np.pi, 0])
	
	counts_guess = sin_fit(bins_r, *popt)-counts_r
	
	plt.plot(bins_r, sin_fit(bins_r, *popt), color='r')
	plt.xlim(left=0,right=1)

	model, params = gen_model(exp_peaks[index], bins_r[peaks], counts_guess[peaks])
	if f == "7mu_Fe.TKA":	# model has trouble fitting, guidance is given
		params["p3_center"].min = 0.229
		params["p3_center"].max = 0.231
		params["p3_amplitude"].max = 20e3
		params["p9_center"].min = 0.743
		params["p9_center"].max = 0.745
		params["p9_amplitude"].max = 20e3
	result = model.fit(counts_guess, params, x=bins_r)
	#plt.plot(bins_r, result.best_fit)
	result_params = list(result.best_values.keys())
	total_fit = np.zeros_like(bins)
	for i in range(exp_peaks[index]):
		amp = result.best_values[result_params[int(i*3)]]
		mu = result.best_values[result_params[int(i*3+1)]]
		sigma = result.best_values[result_params[int(i*3+2)]]
		total_fit += gauss_fit(bins, amp, mu, sigma)
		
	plt.subplot(3, 1, 3)
	plt.stairs(counts_guess, np.arange(n+1)*14620/n*70e-6, color='r', alpha=0.5)
	plt.plot(bins, total_fit)
	plt.xlim(left=0,right=1)
plt.show()
