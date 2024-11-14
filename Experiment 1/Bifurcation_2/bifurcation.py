import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
import os

def bifurcate():
	plt.figure(0)
	sz = 5
	plt.text(0.2,3.25,"Period 1")
	plt.scatter(np.arange(100,3000,100)/1e3,period1,color='k',marker='.',s=sz)
	plt.plot([3,3],[0,3.5],color='k',linestyle='--')
	
	plt.text(4.7,3.25,"Period 2")
	for i in range(2):
		plt.scatter(np.arange(3000,7800,100)/1e3,period2[i],color='k',marker='.',s=sz)
	plt.plot([7.8,7.8],[0,3.5],color='k',linestyle='--')
	
	plt.text(7.9,3.25,"P4")
	for i in range(4):
		plt.scatter(np.arange(7800,8500,100)/1e3,period4[i],color='k',marker='.',s=sz)
	plt.scatter(np.ones(8)*8.5,period8,color='k',marker='.',s=sz)
	plt.plot([8.5,8.5],[0,3.5],color='k',linestyle='--')
	
	plt.text(9.2,3.25,"Chaos 1")
	for volt,chao in zip(np.arange(8600,11600,100)/1e3,chaos):
		plt.scatter(np.ones(len(chao))*volt,list(chao),color='k',marker='.',s=sz)
	plt.plot([11.6,11.6],[0,3.5],color='k',linestyle='--')
	
	plt.text(12.2,3.25,"Period 3")
	for i in range(3):
		plt.scatter(np.arange(11600,14300,100)/1e3,period3[i],color='k',marker='.',s=sz)
	plt.plot([14.3,14.3],[0,3.5],color='k',linestyle='--')
	plt.text(14.6,3.25,"P6")
	for i in range(6):
		plt.scatter(np.arange(14300,15400,100)/1e3,period6[i],color='k',marker='.',s=sz)
	plt.scatter(np.ones(12)*15.4,period12,color='k',marker='.',s=sz)
	plt.plot([15.4,15.4],[0,3.5],color='k',linestyle='--')
	
	plt.text(16.1,3.25,"Chaos 2")
	for volt,chao in zip(np.arange(15500,20100,100)/1e3,chaos2):
		plt.scatter(np.ones(len(chao))*volt,list(chao),color='k',marker='.',s=sz)
	plt.xlabel("$V_{dr}$ ($\\mathrm{V_{PP}}$)")
	plt.ylabel("$V_{R}$ (V)")
	plt.show()

def quickplot(start,end):
	for i in range(start,end,100):
		peaks = []
		for j in range(1,33):
			time,Vdrive,VI = np.loadtxt(f'./{i}/{i}_{j:02}.csv',delimiter=',',skiprows=3).T
			
			VI_smooth = savgol_filter(VI,51,3)
			guess_peaks,_ = find_peaks(VI_smooth*100,distance=500)
			
			peaks += [max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks]
		if i <= 3000:
			peaks = [peak/1000 for peak in peaks]
		plt.scatter(np.ones(len(peaks))*i,peaks,color='k',marker='.')
	#plt.plot(np.arange(start,end), np.arange(start,end)/1000*0.15-0.1)
	#plt.plot(np.arange(start,end), np.arange(start,end)/1000*0.13-0.45)
	plt.show()

period1 = []
period2 = []
period4 = []
period8 = []
chaos = []
period3 = []
period6 = []
period12 = []
chaos2 = []

if not os.path.isfile("period1.npy"):
	for i in range(100,3000,100):
		peaks = []
		for j in range(1,33):
			# extract columns from CSV, convert to row arrays
			time,Vdrive,VI = np.loadtxt(f'./{i}/{i}_{j:02}.csv',delimiter=',',skiprows=3).T
			
			# smooth out the VI curve to find guess the location of peaks
			VI_smooth = savgol_filter(VI,51,3)
			guess_peaks,_ = find_peaks(VI_smooth*100,distance=500)
			
			# obtain actual values of peak, append to existing list
			peaks += [max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks]
			
		mean_peak = np.mean(peaks)
		# data for Vdrive = 3 V are in mV, so need to scale values
		if i <= 3000:
			mean_peak = mean_peak/1000
		period1.append(mean_peak)
		# array of x = Vdrive for making scatterplot
	period1 = np.array(period1)
	np.save("period1",period1)
else:
	period1 = np.load("period1.npy")

if not os.path.isfile("period2.npy"):
	for i in range(3000,7800,100):
		peaks = []
		for j in range(1,33):
			time,Vdrive,VI = np.loadtxt(f'./{i}/{i}_{j:02}.csv',delimiter=',',skiprows=3).T
			VI_smooth = savgol_filter(VI,51,3)
			guess_peaks,_ = find_peaks(VI_smooth*100,distance=500)
			peaks += [max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks]
		
		midpoint = np.mean(peaks)
		
		mean_peak_2 = np.mean([peak for peak in peaks if peak > midpoint])
		mean_peak_1 = np.mean([peak for peak in peaks if peak < midpoint])
		# data for Vdrive = 3 V are in mV, so need to scale values
		if i <= 3000:
			mean_peak_2 = mean_peak_2/1000
			mean_peak_1 = mean_peak_1/1000
		period2.append([mean_peak_1,mean_peak_2])
	period2 = np.array(period2).T
	np.save("period2",period2)
else:
	period2 = np.load("period2.npy")

if not os.path.isfile("period4.npy"):
	for i in range(7800,8500,100):
		peaks = []
		for j in range(1,33):
			time,Vdrive,VI = np.loadtxt(f'./{i}/{i}_{j:02}.csv',delimiter=',',skiprows=3).T
			VI_smooth = savgol_filter(VI,51,3)
			guess_peaks,_ = find_peaks(VI_smooth*100,distance=500)
			peaks += [max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks]

		midpoint = np.mean(peaks)
		upper_quarter = np.mean([peak for peak in peaks if peak > midpoint])
		lower_quarter = np.mean([peak for peak in peaks if peak < midpoint])
		
		mean_peak_4 = np.mean([peak for peak in peaks if peak > upper_quarter])
		mean_peak_3 = np.mean([peak for peak in peaks if upper_quarter > peak > midpoint])
		mean_peak_2 = np.mean([peak for peak in peaks if midpoint > peak > lower_quarter])
		mean_peak_1 = np.mean([peak for peak in peaks if lower_quarter > peak])
		period4.append([mean_peak_1,mean_peak_2,mean_peak_3,mean_peak_4])
	period4 = np.array(period4).T
	np.save("period4",period4)
else:
	period4 = np.load("period4.npy")

if not os.path.isfile("period8.npy"):
	i = 8500
	peaks = []
	for j in range(1,33):
		time,Vdrive,VI = np.loadtxt(f'./{i}/{i}_{j:02}.csv',delimiter=',',skiprows=3).T
		
		VI_smooth = savgol_filter(VI,51,3)
		guess_peaks,_ = find_peaks(VI_smooth*100,distance=500)
		
		peaks += [max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks]
		
	unique_peaks,counts = np.unique(peaks,return_counts=True,equal_nan=False)
	period8 = [ np.average(unique_peaks[0:2],weights=counts[0:2]),
				np.average(unique_peaks[2:4],weights=counts[2:4]),
				np.average(unique_peaks[4:6],weights=counts[4:6]),
				unique_peaks[6],
				np.average(unique_peaks[7:10],weights=counts[7:10]),
				np.average(unique_peaks[10:12],weights=counts[10:12]),
				unique_peaks[12],
				np.average(unique_peaks[13:],weights=counts[13:])]
	period8 = np.array(period8)
	np.save("period8",period8)
else:
	period8 = np.load("period8.npy")

if not os.path.isfile("chaos.npy"):
	for i in range(8600,11600,100):
		peaks = set()
		for j in range(1,33):
			time,Vdrive,VI = np.loadtxt(f'./{i}/{i}_{j:02}.csv',delimiter=',',skiprows=3).T
			
			VI_smooth = savgol_filter(VI,51,3)
			guess_peaks,_ = find_peaks(VI_smooth*100,distance=500)
			
			peaks.update([max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks])
		chaos.append(peaks)
	np.save("chaos",chaos)
else:
	chaos = np.load("chaos.npy",allow_pickle=True)

if not os.path.isfile("period3.npy"):
	for i in range(11600,14300,100):
		peaks = []
		for j in range(1,33):
			time,Vdrive,VI = np.loadtxt(f'./{i}/{i}_{j:02}.csv',delimiter=',',skiprows=3).T
			
			VI_smooth = savgol_filter(VI,51,3)
			guess_peaks,_ = find_peaks(VI_smooth*100,distance=500)
			
			peaks += [max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks]
		upper_third = i/1000*0.15-0.1
		lower_third = i/1000*0.13-0.45
		mean_peak_3 = np.mean([peak for peak in peaks if peak > upper_third])
		mean_peak_2 = np.mean([peak for peak in peaks if upper_third > peak > lower_third])
		mean_peak_1 = np.mean([peak for peak in peaks if lower_third > peak])
		period3.append([mean_peak_1,mean_peak_2,mean_peak_3])
	period3 = np.array(period3).T
	np.save("period3",period3)
else:
	period3 = np.load("period3.npy")

if not os.path.isfile("period6.npy"):
	for i in range(14300,15400,100):
		peaks = []
		for j in range(1,33):
			time,Vdrive,VI = np.loadtxt(f'./{i}/{i}_{j:02}.csv',delimiter=',',skiprows=3).T
			
			VI_smooth = savgol_filter(VI,51,3)
			guess_peaks,_ = find_peaks(VI_smooth*100,distance=500)
			
			peaks += [max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks]
		upper_third = i/1000*0.15-0.1
		lower_third = i/1000*0.13-0.45
		upper_sixth = np.mean([peak for peak in peaks if peak > upper_third])
		inter_sixth = np.mean([peak for peak in peaks if upper_third > peak > lower_third])
		lower_sixth = np.mean([peak for peak in peaks if lower_third > peak])
		
		mean_peak_6 = np.mean([peak for peak in peaks if peak > upper_sixth])
		mean_peak_5 = np.mean([peak for peak in peaks if upper_sixth > peak > upper_third])
		mean_peak_4 = np.mean([peak for peak in peaks if upper_third > peak > inter_sixth])
		mean_peak_3 = np.mean([peak for peak in peaks if inter_sixth > peak > lower_third])
		mean_peak_2 = np.mean([peak for peak in peaks if lower_third > peak > lower_sixth])
		mean_peak_1 = np.mean([peak for peak in peaks if lower_sixth > peak])
		period6.append([mean_peak_1,mean_peak_2,mean_peak_3,mean_peak_4,mean_peak_5,mean_peak_6])
	period6 = np.array(period6).T
	np.save("period6",period6)
else:
	period6 = np.load("period6.npy")
if not os.path.isfile("period12.npy"):
	i = 15400	
	peaks = []
	for j in range(1,33):
		time,Vdrive,VI = np.loadtxt(f'./{i}/{i}_{j:02}.csv',delimiter=',',skiprows=3).T
		
		VI_smooth = savgol_filter(VI,51,3)
		guess_peaks,_ = find_peaks(VI_smooth*100,distance=500)
		
		peaks += [max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks]
	unique_peaks,counts = np.unique(peaks,return_counts=True,equal_nan=False)
	period12 = [unique_peaks[0],unique_peaks[1],unique_peaks[2],unique_peaks[3],
				np.average(unique_peaks[4:6],weights=counts[4:6]),
				np.average(unique_peaks[6:8],weights=counts[6:8]),
				np.average(unique_peaks[8:10],weights=counts[8:10]),
				unique_peaks[10],
				np.average(unique_peaks[11:13],weights=counts[11:13]),
				np.average(unique_peaks[13:15],weights=counts[13:15]),
				unique_peaks[15], unique_peaks[16]]
	period12 = np.array(period12)
	np.save("period12",period12)
else:
	period12 = np.load("period12.npy")
	
if not os.path.isfile("chaos2.npy"):
	for i in range(15500,20100,100):
		peaks = set()
		for j in range(1,33):
			time,Vdrive,VI = np.loadtxt(f'./{i}/{i}_{j:02}.csv',delimiter=',',skiprows=3).T
			
			VI_smooth = savgol_filter(VI,51,3)
			guess_peaks,_ = find_peaks(VI_smooth*100,distance=500)
			
			peaks.update([max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks])
		chaos2.append(peaks)
	np.save("chaos2",chaos2)
else:
	chaos2 = np.load("chaos2.npy",allow_pickle=True)

bifurcate()
