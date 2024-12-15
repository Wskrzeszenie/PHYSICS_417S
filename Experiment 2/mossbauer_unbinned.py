from mossbauer_common import *

file_check()

for index, f in enumerate(filelist):
	plt.figure(index)
	plt.suptitle(f)
	counts = np.loadtxt(f+".TKA")[2:14262]	# so that data is 14260 = 2*2*5*23*31, can rebin into nice bin sizes
	bins = np.arange(len(counts))*70e-6+35e-6
	
	plt.subplot(3, 1, 1)
	plt.stairs(counts, np.arange(len(counts)+1)*70e-6, color='k')
	plt.xlim(left=0,right=1)
	
	smooth = sp.signal.savgol_filter(-counts, 51, 3)
	prom_guess = 0
	peaks, _ = sp.signal.find_peaks(smooth, distance=100, prominence=prom_guess)
	while len(peaks) != exp_peaks[index]:
		prom_guess += 1
		peaks, _ = sp.signal.find_peaks(smooth, distance=100, prominence=prom_guess)
		
	counts_filtered = np.copy(counts)
	bins_filtered = np.copy(bins)
	w = 100
	print(f, len(peaks))
	for peak in peaks:
		counts_filtered[(peak-w):(peak+w)] = np.nan
		bins_filtered[(peak-w):(peak+w)] = np.nan
	
	plt.subplot(3, 1, 2)
	plt.stairs(counts_filtered, np.arange(len(counts)+1)*70e-6, color='k')
	
	counts_filtered = counts_filtered[~np.isnan(counts_filtered)]
	bins_filtered = bins_filtered[~np.isnan(bins_filtered)]
	
	popt, pcov = opt.curve_fit(exact_fit, bins_filtered, counts_filtered, 
								p0=[np.mean(counts_filtered), np.mean(counts_filtered), (max(counts_filtered)-min(counts_filtered))/2, 0])
	print(popt[-1])
	popt, pcov = opt.curve_fit(sin_fit, bins_filtered, counts_filtered, 
								p0=[-np.mean(counts_filtered), (max(counts_filtered)-min(counts_filtered))/2, 0])
	print(popt[-1])
	popt, pcov = opt.curve_fit(sin_fit, bins_filtered, counts_filtered, 
								p0=[np.mean(counts_filtered), (max(counts_filtered)-min(counts_filtered))/2, 0])
	print(popt[-1])
	perr = np.sqrt(np.diag(pcov))
	
	plt.plot(bins, sin_fit(bins, *popt), color='r')
	plt.xlim(left=0,right=1)
	
	counts_cleaned = sin_fit(bins, *popt)-counts
	for peak in peaks:
		print(int(max(counts_cleaned[(peak-w):(peak+w)])), end=' ')
	print()
	
	model, params = gen_model(exp_peaks[index], bins[peaks], counts_cleaned[peaks])
	result = model.fit(counts_cleaned, params, x=bins)
	vE_data = calculate_energies(result, popt, perr, index)
	np.save(f+"_vE_data", vE_data)
	plt.subplot(3, 1, 3)
	plt.stairs(counts_cleaned/max(counts_cleaned), np.arange(len(counts)+1)*70e-6)
	for peak in peaks:
		plt.plot([bins[peak-w],bins[peak-w]], [0,1], color='k', linestyle='--')
		plt.plot([bins[peak+w],bins[peak+w]], [0,1], color='k', linestyle='--')
	plt.plot(bins, result.best_fit/max(counts_cleaned))
	plt.xlim(left=0,right=1)
	
plt.show()