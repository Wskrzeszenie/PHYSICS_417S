import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
import os

np.set_printoptions(linewidth=np.inf)

def bifurcate():
	plt.figure(0)
	sz = 10
	for i in range(2):
		plt.scatter(np.arange(7800,7830,10)/1e3,period2[i],color='k',marker='.',s=sz)
	for i in range(4):
		plt.scatter(np.arange(7830,8460,10)/1e3,period4[i],color='k',marker='.',s=sz)
	for i in range(8):
		plt.scatter(np.arange(8460,8580,10)/1e3,period8[i],color='k',marker='.',s=sz)
	for volt,chao in zip(np.arange(8580,8690,10)/1e3,chaos):
		plt.scatter(np.ones(len(chao))*volt,list(chao),color='k',marker='.',s=sz)
	plt.scatter(np.ones(12)*8.69,period3_early[0],color='k',marker='.',s=sz)
	plt.scatter(np.ones(12)*8.7,period3_early[1],color='k',marker='.',s=sz)
	for volt,chao in zip(np.arange(8710,9010,10)/1e3,chaos2):
		plt.scatter(np.ones(len(chao))*volt,list(chao),color='k',marker='.',s=sz)
	plt.xlabel("$V_{dr}$ ($\\mathrm{V_{PP}}$)")
	plt.ylabel("$V_{R}$ (V)")
	
def bifurcate2():
	plt.figure(1)
	sz = 10
	for i in range(3):
		plt.scatter(np.arange(14300,14420,10)/1e3,period3[i],color='k',marker='.',s=sz)
	for i in range(6):
		plt.scatter(np.arange(14420,15350,10)/1e3,period6[i],color='k',marker='.',s=sz)
	for i in range(12):
		plt.scatter(np.arange(15350,15760,10)/1e3,period12[i],color='k',marker='.',s=sz)
	for volt,chao in zip(np.arange(15760,16010,10)/1e3,chaos3):
		plt.scatter(np.ones(len(chao))*volt,list(chao),color='k',marker='.',s=sz)
	plt.xlabel("$V_{dr}$ ($\\mathrm{V_{PP}}$)")
	plt.ylabel("$V_{R}$ (V)")

def split3(x):
	if x < 14500:
		return 2.46
	elif 14500 <= x <= 15900:
		return (2.66-2.46)/(15900-14500)*(x-14500) + 2.46
	return 2.66

def split2(x):
	if x < 14500:
		return 1.76
	elif 14500 <= x <= 15990:
		return (1.99-1.76)/(15990-14500)*(x-14500) + 1.76
	return 1.99

def split1(x):
	if x < 14520:
		return 1.2
	elif 14520 <= x < 14830:
		return 1.24
	elif 14830 <= x < 15230:
		return 1.28
	elif 15230 <= x < 15530:
		return 1.32
	elif 15530 <= x < 15810:
		return 1.36
	return 1.4

def quickplot(start,end):
	for i in range(start,end,10):
		peaks = []
		for j in range(1,33):
			time,Vdrive,VI = np.loadtxt(f'./{i}/{i}_{j:02}.csv',delimiter=',',skiprows=3).T
			
			VI_smooth = savgol_filter(VI,51,3)
			guess_peaks,_ = find_peaks(VI_smooth,distance=500)
			
			peaks += [max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks]
		plt.scatter(np.ones(len(peaks))*i,peaks,color='k',marker='.')
	plt.show()
	
def quickplot_chaos(start,end):
	for i in range(start,end,10):
		peaks = []
		for j in range(1,2):
			time,Vdrive,VI = np.loadtxt(f'./{i}/{i}_{j:02}.csv',delimiter=',',skiprows=3).T
			
			VI_smooth = savgol_filter(VI,51,3)
			guess_peaks,_ = find_peaks(VI_smooth,distance=500)
			
			peaks += [max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks]
		plt.scatter(np.ones(len(peaks))*i,peaks,color='k',marker='.')
	plt.plot([start,end], [2.25,2.25])
	plt.plot([start,end], [1.6,1.6])
	plt.plot(np.arange(start,end,10),[split3(x) for x in np.arange(start,end,10)])
	plt.plot(np.arange(start,end,10),[split2(x) for x in np.arange(start,end,10)])
	plt.plot(np.arange(start,end,10),[split1(x) for x in np.arange(start,end,10)])
	plt.show()

period2 = []
period4 = []
period8 = []
chaos = []
period3_early = []
chaos2 = []
period3 = []
period6 = []
period12 = []
chaos3 = []

#quickplot(7800,9010)
#quickplot(14300,16510)

if not os.path.isfile("period2.npy"):
	for i in range(7800,7830,10):
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
		period2.append([mean_peak_1,mean_peak_2])
	period2 = np.array(period2).T
	np.save("period2",period2)
else:
	period2 = np.load("period2.npy")

if not os.path.isfile("period4.npy"):
	for i in range(7830,8460,10):
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
	for i in range(8460,8580,10):
		peaks = []
		for j in range(1,33):
			time,Vdrive,VI = np.loadtxt(f'./{i}/{i}_{j:02}.csv',delimiter=',',skiprows=3).T
			VI_smooth = savgol_filter(VI,51,3)
			guess_peaks,_ = find_peaks(VI_smooth*100,distance=500)
			peaks += [max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks]

		midpoint = np.mean(peaks)
		upper_quarter = np.mean([peak for peak in peaks if peak > midpoint])
		lower_quarter = np.mean([peak for peak in peaks if peak < midpoint])
		
		uu_eighth = np.mean([peak for peak in peaks if peak > upper_quarter])
		lu_eighth = np.mean([peak for peak in peaks if upper_quarter > peak > midpoint])
		ul_eighth = np.mean([peak for peak in peaks if midpoint > peak > lower_quarter])
		ll_eighth = np.mean([peak for peak in peaks if lower_quarter > peak])
		
		mean_peak_8 = np.mean([peak for peak in peaks if peak > uu_eighth])
		mean_peak_7 = np.mean([peak for peak in peaks if uu_eighth > peak > upper_quarter])
		mean_peak_6 = np.mean([peak for peak in peaks if upper_quarter > peak > lu_eighth])
		mean_peak_5 = np.mean([peak for peak in peaks if lu_eighth > peak > midpoint])
		mean_peak_4 = np.mean([peak for peak in peaks if midpoint > peak > ul_eighth])
		mean_peak_3 = np.mean([peak for peak in peaks if ul_eighth > peak > lower_quarter])
		mean_peak_2 = np.mean([peak for peak in peaks if lower_quarter > peak > ll_eighth])
		mean_peak_1 = np.mean([peak for peak in peaks if ll_eighth > peak])
		
		period8.append([mean_peak_1,mean_peak_2,mean_peak_3,mean_peak_4,mean_peak_5,mean_peak_6,mean_peak_7,mean_peak_8])
	period8 = np.array(period8).T
	np.save("period8",period8)
else:
	period8 = np.load("period8.npy")
if not os.path.isfile("chaos.npy"):
	for i in range(8580,8690,10):
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

if not os.path.isfile("period3_early.npy"):
	i = 8690
	peaks = []
	for j in range(1,33):
		time,Vdrive,VI = np.loadtxt(f'./{i}/{i}_{j:02}.csv',delimiter=',',skiprows=3).T
		VI_smooth = savgol_filter(VI,51,3)
		guess_peaks,_ = find_peaks(VI_smooth*100,distance=500)
		peaks += [max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks]
		
	unique_peaks,counts = np.unique(peaks,return_counts=True,equal_nan=False)
	period3_early.append([  unique_peaks[0],
							np.average(unique_peaks[1:3],weights=counts[1:3]),
							np.average(unique_peaks[3:6],weights=counts[3:6]),
							np.average(unique_peaks[6:8],weights=counts[6:8]),
							unique_peaks[8],
							np.average(unique_peaks[9:11],weights=counts[9:11]),
							np.average(unique_peaks[11:13],weights=counts[11:13]),
							np.average(unique_peaks[13:15],weights=counts[13:15]),
							np.average(unique_peaks[15:19],weights=counts[15:19]),
							np.average(unique_peaks[19:22],weights=counts[19:22]),
							unique_peaks[22],unique_peaks[23]])
	i = 8700
	peaks = []
	for j in range(1,33):
		time,Vdrive,VI = np.loadtxt(f'./{i}/{i}_{j:02}.csv',delimiter=',',skiprows=3).T
		VI_smooth = savgol_filter(VI,51,3)
		guess_peaks,_ = find_peaks(VI_smooth*100,distance=500)
		peaks += [max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks]
		
	unique_peaks,counts = np.unique(peaks,return_counts=True,equal_nan=False)
	#print(unique_peaks[unique_peaks < 1.0])
	#print(counts[np.where(unique_peaks < 1.0)])
	period3_early.append([  unique_peaks[0],
							np.average(unique_peaks[1:3],weights=counts[1:3]),
							unique_peaks[3],
							np.average(unique_peaks[4:6],weights=counts[4:6]),
							unique_peaks[6],unique_peaks[7],
							np.average(unique_peaks[8:10],weights=counts[8:10]),
							np.average(unique_peaks[10:12],weights=counts[10:12]),
							unique_peaks[12],
							np.average(unique_peaks[13:15],weights=counts[13:15]),
							np.average(unique_peaks[15:17],weights=counts[15:17]),
							np.average(unique_peaks[17:19],weights=counts[17:19])])
	np.save("period3_early",period3_early)
else:
	period3_early = np.load("period3_early.npy")
#np.save("chaos",chaos)

if not os.path.isfile("chaos2.npy"):
	for i in range(8710,9010,10):
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

#bifurcate()

if not os.path.isfile("period3.npy"):
	for i in range(14300,14420,10):
		peaks = []
		for j in range(1,33):
			time,Vdrive,VI = np.loadtxt(f'./{i}/{i}_{j:02}.csv',delimiter=',',skiprows=3).T
			VI_smooth = savgol_filter(VI,51,3)
			guess_peaks,_ = find_peaks(VI_smooth*100,distance=500)
			peaks += [max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks]
		upper_third = 2.25
		lower_third = 1.6
		
		mean_peak_3 = np.mean([peak for peak in peaks if upper_third < peak])
		mean_peak_2 = np.mean([peak for peak in peaks if lower_third < peak < upper_third])
		mean_peak_1 = np.mean([peak for peak in peaks if peak < lower_third])
		period3.append([mean_peak_1,mean_peak_2,mean_peak_3])
	period3 = np.array(period3).T
	np.save("period3",period3)
else:
	period3 = np.load("period3.npy")

if not os.path.isfile("period6.npy"):
	for i in range(14420,15350,10):
		peaks = []
		for j in range(1,33):
			time,Vdrive,VI = np.loadtxt(f'./{i}/{i}_{j:02}.csv',delimiter=',',skiprows=3).T
			VI_smooth = savgol_filter(VI,51,3)
			guess_peaks,_ = find_peaks(VI_smooth*100,distance=500)
			peaks += [max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks]
		upper_third = 2.25
		lower_third = 1.6
		
		mean_peak_6 = np.mean([peak for peak in peaks if split3(i) < peak])
		mean_peak_5 = np.mean([peak for peak in peaks if upper_third < peak < split3(i)])
		mean_peak_4 = np.mean([peak for peak in peaks if split2(i) < peak < upper_third])
		mean_peak_3 = np.mean([peak for peak in peaks if lower_third < peak < split2(i)])
		mean_peak_2 = np.mean([peak for peak in peaks if split1(i) < peak < lower_third])
		mean_peak_1 = np.mean([peak for peak in peaks if peak < split1(i)])
		period6.append([mean_peak_1,mean_peak_2,mean_peak_3,mean_peak_4,mean_peak_5,mean_peak_6])
	period6 = np.array(period6).T
	np.save("period6",period6)
else:
	period6 = np.load("period6.npy")

if not os.path.isfile("period12.npy"):
	for i in range(15350,15760,10):
		peaks = []
		for j in range(1,33):
			time,Vdrive,VI = np.loadtxt(f'./{i}/{i}_{j:02}.csv',delimiter=',',skiprows=3).T
			VI_smooth = savgol_filter(VI,51,3)
			guess_peaks,_ = find_peaks(VI_smooth*100,distance=500)
			peaks += [max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks]
		upper_third = 2.25
		lower_third = 1.6

		u3_12 = np.mean([peak for peak in peaks if split3(i) < peak])
		u2_12 = np.mean([peak for peak in peaks if upper_third < peak < split3(i)])
		u1_12 = np.mean([peak for peak in peaks if split2(i) < peak < upper_third])
		l3_12 = np.mean([peak for peak in peaks if lower_third < peak < split2(i)])
		l2_12 = np.mean([peak for peak in peaks if split1(i) < peak < lower_third])
		l1_12 = np.mean([peak for peak in peaks if peak < split1(i)])
		
		# fix irregularity of spliting due to discretization of values
		if i == 15530: l3_12 = 1.79
		if i == 15640: l2_12 = 1.44
		if i == 15660: u1_12 = 2.03
		
		mean_peak_12 = np.mean([peak for peak in peaks if u3_12 < peak])
		mean_peak_11 = np.mean([peak for peak in peaks if split3(i) < peak < u3_12])
		mean_peak_10 = np.mean([peak for peak in peaks if u2_12 < peak < split3(i)])
		mean_peak_9 = np.mean([peak for peak in peaks if upper_third < peak < u2_12])
		mean_peak_8 = np.mean([peak for peak in peaks if u1_12 < peak < upper_third])
		mean_peak_7 = np.mean([peak for peak in peaks if split2(i) < peak < u1_12])
		mean_peak_6 = np.mean([peak for peak in peaks if l3_12 < peak < split2(i)])
		mean_peak_5 = np.mean([peak for peak in peaks if lower_third < peak < l3_12])
		mean_peak_4 = np.mean([peak for peak in peaks if l2_12 < peak < lower_third])
		mean_peak_3 = np.mean([peak for peak in peaks if split1(i) < peak < l2_12])
		mean_peak_2 = np.mean([peak for peak in peaks if l1_12 < peak < split1(i)])
		mean_peak_1 = np.mean([peak for peak in peaks if peak < l1_12])
		
		period12.append([mean_peak_1,mean_peak_2,mean_peak_3,mean_peak_4,mean_peak_5,mean_peak_6,mean_peak_7,mean_peak_8,mean_peak_9,mean_peak_10,mean_peak_11,mean_peak_12])
	period12 = np.array(period12).T
	np.save("period12",period12)
else:
	period12 = np.load("period12.npy")

if not os.path.isfile("chaos3.npy"):
	for i in range(15760,16010,10):
		peaks = set()
		for j in range(1,33):
			time,Vdrive,VI = np.loadtxt(f'./{i}/{i}_{j:02}.csv',delimiter=',',skiprows=3).T
			VI_smooth = savgol_filter(VI,51,3)
			guess_peaks,_ = find_peaks(VI_smooth*100,distance=500)
			peaks.update([max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks])
		chaos3.append(peaks)
	np.save("chaos3",chaos3)
else:
	chaos3 = np.load("chaos3.npy",allow_pickle=True)

bifurcate()
bifurcate2()
plt.show()