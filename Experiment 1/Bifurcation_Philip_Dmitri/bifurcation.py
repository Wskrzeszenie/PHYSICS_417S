import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.signal import find_peaks

plt.figure(0)

for i in range(100,20100,100):
	
	# extract columns from CSV, convert to row arrays
	time, Vdrive, VI = np.loadtxt((f'{i:05}' + "V.csv"),delimiter=',',skiprows=3).T
	
	# smooth out the VI curve to find guess the location of peaks
	VI_smooth = savgol_filter(VI,51,3)
	guess_peaks, _ = find_peaks(VI_smooth*100, distance=500)
	
	# obtain actual values of peak
	peaks = np.array([max(VI[max(0,i-50):min(len(VI)-1,i+50)]) for i in guess_peaks])
	
	# data for Vdrive = 3 V are in mV, so need to scale values
	if i <= 3000:
		peaks = peaks/1000
	# array of x = Vdrive for making scatterplot
	X = [i/1000 for _ in range(len(peaks))]
	plt.scatter(X,peaks,color='k',marker='.')
#plt.xlim(7.5,10)
#plt.ylim(0.5,2)
plt.show()
