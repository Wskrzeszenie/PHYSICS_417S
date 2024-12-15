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
	
	n = 14260	# 1 2 4 5 10 20 23 31 46 62 92 115 124 155 230 310 460 620 713 1426 2852 3565 7130 14260
	bins_r = np.arange(n)*14620/n*70e-6+14620/n*35e-6
	counts_r = np.zeros(n)
	for i in range(14260):
		counts_r[i//(14260//n)] += counts[i]
	
	popt, pcov = opt.curve_fit(sin_fit, bins_r, counts_r, p0=[np.mean(counts_r), (max(counts_r)-min(counts_r))/2, 0])
	
	counts_guess = sin_fit(bins_r, *popt)-counts_r
	prom_guess = 0
	peaks, _ = sp.signal.find_peaks(counts_guess/max(counts_guess)*100, distance = 5*n/460, prominence=prom_guess)
	while len(peaks) != exp_peaks[index]:
		prom_guess += 1
		peaks, _ = sp.signal.find_peaks(counts_guess/max(counts_guess)*100, distance = 5*n/460, prominence=prom_guess)

	counts_rf = np.copy(counts_r)
	bins_rf = np.copy(bins_r)
	w = int(5*n/460)
	print(f, len(peaks))
	for peak in peaks:
		counts_rf[(peak-w):(peak+w)] = np.nan
		bins_rf[(peak-w):(peak+w)] = np.nan
		
	plt.subplot(3, 1, 2)
	plt.stairs(counts_rf, np.arange(n+1)*14620/n*70e-6, color='k')
	
	counts_rf = counts_rf[~np.isnan(counts_rf)]
	bins_rf = bins_rf[~np.isnan(bins_rf)]
	
	popt, pcov = opt.curve_fit(sin_fit, bins_rf, counts_rf, 
								p0=[np.mean(counts_rf), (max(counts_rf)-min(counts_rf))/2, 0])
	perr = np.sqrt(np.diag(pcov))
	
	plt.plot(bins_r, sin_fit(bins_r, *popt), color='r')
	plt.xlim(left=0,right=1)
	
	counts_cleaned = sin_fit(bins_r, *popt)-counts_r
	for peak in peaks:
		print(max(counts_filtered[(peak-w):(peak+w)]))
		
	model, params = gen_model(exp_peaks[index], bins_r[peaks], counts_cleaned[peaks])
	result = model.fit(counts_cleaned, params, x=bins_r)
	vE_data = calculate_energies(result, popt, perr, index)
	np.save("vE_data", vE_data)
	plt.subplot(3, 1, 3)
	plt.stairs(counts_cleaned, np.arange(n+1)*14620/n*70e-6, color='r', alpha=0.5)
	for peak in peaks:
		plt.plot([bins[peak-w],bins[peak-w]], [0,1], color='k', linestyle='--')
		plt.plot([bins[peak+w],bins[peak+w]], [0,1], color='k', linestyle='--')
	plt.plot(bins_r, result.best_fit)
	plt.xlim(left=0,right=1)
	
#plt.show()
