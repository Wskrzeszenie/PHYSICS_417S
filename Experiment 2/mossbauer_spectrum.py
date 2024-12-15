from mossbauer_common import *

file_check()

for index, f in enumerate(filelist):
	plt.figure(index)
	counts = np.loadtxt(f+".TKA")[2:14262]	# so that data is 14260 = 2*2*5*23*31, can rebin into nice bin sizes
	bins = np.arange(len(counts))*70e-6+35e-6
	
	smooth = sp.signal.savgol_filter(-counts, 51, 3)
	prom_guess = 0
	peaks, _ = sp.signal.find_peaks(smooth, distance=100, prominence=prom_guess)
	while len(peaks) != exp_peaks[index]:
		prom_guess += 1
		peaks, _ = sp.signal.find_peaks(smooth, distance=100, prominence=prom_guess)
		
	counts_filtered = np.copy(counts)
	bins_filtered = np.copy(bins)
	w = 100
	
	for peak in peaks:
		counts_filtered[(peak-w):(peak+w)] = np.nan
		bins_filtered[(peak-w):(peak+w)] = np.nan
	
	counts_filtered = counts_filtered[~np.isnan(counts_filtered)]
	bins_filtered = bins_filtered[~np.isnan(bins_filtered)]
	
	popt, pcov = opt.curve_fit(sin_fit, bins_filtered, counts_filtered, 
								p0=[np.mean(counts_filtered), (max(counts_filtered)-min(counts_filtered))/2, 0])
	
	counts_cleaned = counts-sin_fit(bins, *popt)
	v_bins = v(bins, popt[-1])
	#plt.scatter(v_bins, counts_cleaned/max(counts_cleaned), s=0.1)
	plt.stairs(-counts_cleaned, np.arange(len(counts)+1)*70e-6, color='k')
	plt.title(latex_names[index])
	plt.xlabel("t (s)")
	plt.ylabel("Counts")
	plt.savefig("spectrum_"+f+".png", dpi=300)
plt.show()