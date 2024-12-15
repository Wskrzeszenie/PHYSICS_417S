from mossbauer_common import *

file_check()

# first plot
plt.figure(0)
counts = np.loadtxt("alpha_Fe.TKA")	# so that data is 14260 = 2*2*5*23*31, can rebin into nice bin sizes
counts[:2] = 0
n = 14260
bins = np.arange(n)*70e-6+35e-6

plt.stairs(counts, color='k')

plt.title("Raw Data")
plt.xlabel("Channel")
plt.ylabel("counts")
plt.xlim(left=0)
plt.savefig('example0.png', dpi=300)

# second plot
plt.figure(1)
counts = counts[2:(n+2)]
plt.stairs(counts, np.arange(n+1)*70e-6, color='k')

plt.title("Raw Data (0 < t < 1)")
plt.xlabel("t (s)")
plt.ylabel("Counts")
plt.xlim(left=0,right=1)
plt.savefig('example1.png', dpi=300)

# third plot
plt.figure(2)
smooth = sp.signal.savgol_filter(-counts, 51, 3)
prom_guess = 0
peaks, _ = sp.signal.find_peaks(smooth, distance=100, prominence=prom_guess)
while len(peaks) != 12:
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
print(popt)
perr = np.sqrt(np.diag(pcov))

plt.stairs(counts, np.arange(n+1)*70e-6, color='k')
plt.plot(bins, sin_fit(bins, *popt), color='r')

plt.title("Data w/ Curve Fit for Background")
plt.xlabel("t (s)")
plt.ylabel("Counts")
plt.xlim(left=0,right=1)
plt.savefig('example2.png', dpi=300)

# fourth plot
plt.figure(3)
counts_cleaned = sin_fit(bins, *popt)-counts
plt.stairs(counts_cleaned, np.arange(n+1)*70e-6, color='k')

plt.title("Data w/ Background Removed")
plt.xlabel("t (s)")
plt.ylabel("Counts")
plt.xlim(left=0,right=1)
plt.savefig('example3.png', dpi=300)

# fifth plot
plt.figure(4)
model, params = gen_model(12, bins[peaks], counts_cleaned[peaks])
result = model.fit(counts_cleaned, params, x=bins)
plt.stairs(counts_cleaned, np.arange(n+1)*70e-6, color='k')
plt.plot(bins, result.best_fit, color='r')

plt.title("Data w/ Fitted Gaussians")
plt.xlabel("t (s)")
plt.ylabel("Counts")
plt.xlim(left=0,right=1)
plt.savefig('example4.png', dpi=300)

# sixth plot
plt.figure(5)
plt.stairs(counts_cleaned, np.arange(n+1)*70e-6, color='w')
plt.plot(bins, result.best_fit, color='r')

plt.title("Fitted Gaussians")
plt.xlabel("t (s)")
plt.ylabel("Counts")
plt.xlim(left=0,right=1)
plt.savefig('example5.png', dpi=300)

#plt.show()
