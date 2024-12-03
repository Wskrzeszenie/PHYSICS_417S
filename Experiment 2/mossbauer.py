from mossbauer_common import *

file_check()
#calibrate = [90, 180, 55, 65, 135]
#calibrate = [35, 25, 35, 35, 35, 35, 20]
calibrate = [12, 16, 16, 14, 7, 20, 16, 14, 20]

for index, f in enumerate(filelist):
	#plt.subplot(len(filelist),1,index+1)
	plt.suptitle(f)
	plt.figure(index)
	counts = np.loadtxt(f+".TKA")[2:14262]	# so that data is 14260 = 2*2*5*23*31, can rebin into nice bin sizes
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
	w = 5
	print(f, len(peaks))
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
	calculate_energies(result, popt)
		
	plt.subplot(3, 1, 3)
	plt.stairs(counts_guess, np.arange(n+1)*14620/n*70e-6, color='r', alpha=0.5)
	plt.plot(bins_r, result.best_fit)
	plt.xlim(left=0,right=1)
plt.show()
