import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
import os

for i in range(0,8):
	index = 8500 + i*50
	for j in range(1,17):
		time, Vdrive, VI = np.loadtxt(f'./Thickness/{index}/{index}_{j:02}.csv',delimiter=',',skiprows=3).T
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
		'''
		plt.subplot(2,4,i+1)
		plt.scatter(Vdrive, VI, color='k',marker='.')
		#plt.plot(Vdrive[(Vpi-50):(Vpi+50)], VI[(Vpi-50):(Vpi+50)], marker='o', color='r')
		plt.ylim(bottom=1.0)
plt.show()
