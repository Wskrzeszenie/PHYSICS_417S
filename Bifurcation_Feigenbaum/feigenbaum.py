import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
import os

period2 = np.zeros((2,(195-55)//5))
period4 = np.zeros((4,(195-55)//5))

if not os.path.isfile("period2.npy"):
	for i in range(55,195,5):
		peaks = []
		uniques = set()
		for j in range(1,33):
			time,Vdrive,VI = np.loadtxt(f'./{i}_2/{i}_2_{j:02}.csv',delimiter=',',skiprows=3).T
			uniques.update(VI)
			VI_smooth = savgol_filter(VI,51,3)
			guess_peaks,_ = find_peaks(VI_smooth*100,distance=700*85/i)
			peaks += [max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks]
		
		midpoint = np.mean(peaks)
		
		peak_2 = [peak for peak in peaks if peak > midpoint]
		peak_1 = [peak for peak in peaks if peak < midpoint]
		
		uniques = np.array(sorted(list(uniques)))
		err_dig = min(uniques[1:]-uniques[:-1])/2
		
		gap_avg = np.mean(peak_2) - np.mean(peak_1)
		gap_err = err_dig/np.sqrt(len(peaks))
		#gap_err = np.sqrt(np.var(peak_2,ddof=1)+np.var(peak_1,ddof=1))
		if max(peaks) > 10:
			gap_avg /= 1e3
			gap_err /= 1e3
		period2[0,(i-55)//5] = gap_avg
		period2[1,(i-55)//5] = gap_err
	np.save("period2",period2)
else:
	period2 = np.load("period2.npy")


if not os.path.isfile("period4.npy"):
	for i in range(55,195,5):
		peaks = []
		uniques = set()
		for j in range(1,33):
			time,Vdrive,VI = np.loadtxt(f'./{i}_4/{i}_4_{j:02}.csv',delimiter=',',skiprows=3).T
			uniques.update(VI)
			VI_smooth = savgol_filter(VI,51,3)
			guess_peaks,_ = find_peaks(VI_smooth*100,distance=700*85/i)
			peaks += [max(VI[max(0,k-50):min(len(VI)-1,k+50)]) for k in guess_peaks]

		midpoint = np.mean(peaks)
		upper_quarter = np.mean([peak for peak in peaks if peak > midpoint])
		lower_quarter = np.mean([peak for peak in peaks if peak < midpoint])
		
		peak_4 = [peak for peak in peaks if peak > upper_quarter]
		peak_3 = [peak for peak in peaks if upper_quarter > peak > midpoint]
		peak_2 = [peak for peak in peaks if midpoint > peak > lower_quarter]
		peak_1 = [peak for peak in peaks if lower_quarter > peak]
		
		uniques = np.array(sorted(list(uniques)))
		err_dig = min(uniques[1:]-uniques[:-1])/2

		gap2_avg = np.mean(peak_4) - np.mean(peak_3)
		gap2_err = err_dig/np.sqrt(len(peak_4)+len(peak_3))
		#gap2_err = np.sqrt(np.var(peak_4,ddof=1)+np.var(peak_3,ddof=1))	
		gap1_avg = np.mean(peak_2) - np.mean(peak_1)
		gap1_err = err_dig/np.sqrt(len(peak_2)+len(peak_1))
		#gap1_err = np.sqrt(np.var(peak_2,ddof=1)+np.var(peak_1,ddof=1))
		
		period4[0,(i-55)//5] = gap2_avg
		period4[1,(i-55)//5] = gap2_err
		period4[2,(i-55)//5] = gap1_avg
		period4[3,(i-55)//5] = gap1_err
	np.save("period4",period4)
else:
	period4 = np.load("period4.npy")
	
def error_alpha(A,B,dA,dB):
	return np.sqrt((dB*A/B**2)**2
					+(dA/B)**2)

def error_delta(A,B,C,dA,dB,dC):
	return np.sqrt(  (dA/(C-B))**2
					+(dB*(A-C)/(B-C)**2)**2
					+(dC*(A-B)/(B-C)**2)**2)

doubling = np.loadtxt("bifurcation_data.txt",delimiter=',')
doubling_err = np.copy(doubling)
delta_est = np.copy(doubling)
for i in range(len(doubling)):
	for j in range(1,5):
		if doubling[i,j] > 10:
			doubling_err[i,j] = 0.005
		else:
			doubling_err[i,j] = 0.0025
	
V2 = doubling[:,1]
V4 = doubling[:,2]
V8 = doubling[:,3]
V16 = doubling[:,4]

V2_err = doubling_err[:,1]
V4_err = doubling_err[:,2]
V8_err = doubling_err[:,3]
V16_err = doubling_err[:,4]

delta_est[:,1] = (V4-V2)/(V8-V4)
delta_est[:,2] = error_delta(V2,V4,V8,V2_err,V4_err,V8_err)
delta_est[:,3] = (V8-V4)/(V16-V8)
delta_est[:,4] = error_delta(V4,V8,V16,V4_err,V8_err,V16_err)

delta = 4.669201609102990671853203820466201617258185577475768632745651343004134
alpha = 2.502907875095892822283902873218215786381271376727149977336192056779235

plt.figure(0)
plt.subplot(2,3,1)
plt.plot([50,190],[delta, delta],linestyle='--',color='k')
plt.errorbar(np.arange(50,195,5),delta_est[:,1],yerr=delta_est[:,2],marker='o',markersize=2.5,linestyle='',capsize=5,color='r')
plt.errorbar(np.arange(50,195,5),delta_est[:,3],yerr=delta_est[:,4],marker='o',markersize=2.5,linestyle='',capsize=5,color='g')
plt.legend(["Expected $\delta$","1st Estimate","2nd Estimate"])
plt.title("Estimation of 1st Feigenbaum Number")
plt.xlabel("Frequency (kHz)")
plt.ylabel("$\delta$")

plt.subplot(2,3,4)
plt.scatter(np.arange(50,195,5),abs(delta_est[:,1]-delta)/delta*100,s=2.5,color='r')
plt.scatter(np.arange(50,195,5),abs(delta_est[:,3]-delta)/delta*100,s=2.5,color='g')
plt.legend(["1st Estimate","2nd Estimate"])
plt.title("% Error for Estimation of 1st Feigenbaum Number")
plt.xlabel("Frequency (kHz)")
plt.ylabel("% Error")

plt.subplot(2,3,2)
plt.errorbar(np.arange(55,195,5),period2[0],yerr=period2[1],marker='o',markersize=2.5,linestyle='',capsize=5)
plt.title("Gap of 1st Tine")
plt.xlabel("Frequency (kHz)")
plt.ylabel("V")

plt.subplot(2,3,5)
plt.errorbar(np.arange(55,195,5),period4[0],yerr=period4[1],marker='o',markersize=2.5,linestyle='',capsize=5)
plt.errorbar(np.arange(55,195,5),period4[2],yerr=period4[3],marker='o',markersize=2.5,linestyle='',capsize=5)
plt.legend(["Upper Tine Gap","Lower Tine Gap"])
plt.title("Gap of 2nd Tines")
plt.xlabel("Frequency (kHz)")
plt.ylabel("V")

plt.subplot(2,3,3)
plt.plot([55,195],[alpha,alpha],linestyle='--',color='k')
plt.errorbar(np.arange(55,195,5),period2[0]/period4[0],yerr=error_alpha(period2[0],period4[0],period2[1],period4[1]),marker='o',markersize=2.5,linestyle='',capsize=5,color='b')
plt.errorbar(np.arange(55,195,5),period2[0]/period4[2],yerr=error_alpha(period2[0],period4[2],period2[1],period4[3]),marker='o',markersize=2.5,linestyle='',capsize=5,color='r')
plt.legend(["Expected $\\alpha$","Estimate Using Upper Tine","Estimate Using Lower Tine"])
plt.title("Estimation of 2nd Feigenbaum Number")
plt.xlabel("Frequency (kHz)")
plt.ylabel("$\\alpha$")

plt.subplot(2,3,6)
plt.scatter(np.arange(55,195,5),abs(period2[0]/period4[0]-alpha)/alpha*100,s=2.5,color='b')
plt.scatter(np.arange(55,195,5),abs(period2[0]/period4[2]-alpha)/alpha*100,s=2.5,color='r')
plt.legend(["Estimate Using Upper Tine","Estimate Using Lower Tine"])
plt.title("% Error for Estimation of 2nd Feigenbaum Number")
plt.xlabel("Frequency (kHz)")
plt.ylabel("% Error")
plt.show()

