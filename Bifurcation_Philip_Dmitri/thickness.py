import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
import scipy.optimize as opt
import os

plt.figure(0)
'''
for i in range(1, 10):
	index = 8100+i*500
	time, Vdrive, VI = np.loadtxt((f'{index:05}' + "V.csv"),delimiter=',',skiprows=3).T
	plt.subplot(3,3,i)
	plt.scatter(Vdrive, VI, marker='.')
	#plt.xlim(0,20)
	#plt.ylim(0,4)

plt.figure(1)
'''
Vdrives = np.linspace(8.6, 11.7, 32)
thicknesses = np.zeros(32)

# this calculation is very long
# run only once
if not os.path.isfile("./thickness.npy"):
	for i in range(0,32):
		index = 8600 + i*100
		time, Vdrive, VI = np.loadtxt((f'{index:05}' + "V.csv"),delimiter=',',skiprows=3).T
		#if index == 9100:
		#	plt.figure(1)
		#	plt.plot(time, VI)
		#	plt.figure(0)
		'''
		Vdrive_smooth = savgol_filter(Vdrive,63,3)
		VI_smooth = savgol_filter(VI,63,3)
		guess_top = np.argmax(VI_smooth)
		top_index = np.argmax(VI[max(0,guess_top-50):min(guess_top+50,len(VI)-1)])
		VI_peak = (VI[max(0,guess_top-50):min(guess_top+50,len(VI)-1)])[top_index]
		Vdrive_peak = (Vdrive[max(0,guess_top-50):min(guess_top+50,len(VI)-1)])[top_index]
		
		#plt.scatter(Vdrive_smooth, VI_smooth, marker='.')
		temp = Vdrive_smooth[0]
		'''
		Vpi = np.argmax(VI) # Vdrive_peak_index
		VI_peak = np.max(VI)
		thickness = 0 # 0.03937
		search_width = 0.03937
		for j in range(0,len(VI)):
			if VI[j] < 0.5:
				VI[j] = np.nan
				Vdrive[j] = np.nan
		for k in range(0,len(VI)):
			for j in range(Vpi-10,Vpi+10):
				if round(Vdrive[j]-search_width,5) <= round(Vdrive[k],5) <= round(Vdrive[j]+search_width,5):
					thickness = max(thickness, VI[j] - VI[k])
		thicknesses[i] = thickness
		print(index, search_width, thickness)
		#plt.subplot(4,8,i+1)
		#plt.scatter(Vdrive, VI, marker='.')
		#plt.plot(Vdrive[(Vpi-50):(Vpi+50)], VI[(Vpi-50):(Vpi+50)], marker='o', color='r')
		#plt.ylim(bottom=0)
	#plt.show()
	thicknesses = np.array(thicknesses)
	np.save("thickness.npy", thicknesses)

thicknesses = np.load("thickness.npy")

def model(x, A, gamma):
	return A*(x**gamma)
	
Vonset = 8.595
popt, pcov = opt.curve_fit(model, (Vdrives-Vonset), thicknesses)
perr = np.sqrt(np.diag(pcov))
print(popt, perr)
for i in range(0,32):
	index = 8600 + i*100
	time, Vdrive, VI = np.loadtxt((f'{index:05}' + "V.csv"),delimiter=',',skiprows=3).T
	plt.subplot(4,8,i+1)
	plt.scatter(Vdrive, VI, marker='.')
plt.show()
