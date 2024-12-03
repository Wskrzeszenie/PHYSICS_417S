from mossbauer_common import *

file_check()

#calibrate = [90, 180, 55, 65, 135]
calibrate = [136, 89, 46, 66, 108, 58, 94, 69, 76]


for index, f in enumerate(filelist):
	plt.figure(index)
	counts = np.loadtxt(f+".TKA")[2:14262]	# so that data is 14260 = 2*2*5*23*31, can rebin into nice bin sizes
	bins = np.arange(len(counts))*70e-6+35e-6
	
	plt.subplot(3, 1, 1)
	plt.suptitle(f)
	plt.stairs(counts, np.arange(len(counts)+1)*70e-6, color='k')
	plt.xlim(left=0,right=1)
	
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

	model, params = gen_model(exp_peaks[index], bins[peaks], counts_cleaned[peaks])
	result = model.fit(counts_cleaned, params, x=bins)
	calculate_energies(result, popt)
	
	plt.plot(bins, result.best_fit)
	plt.xlim(left=0,right=1)
	
plt.show()
