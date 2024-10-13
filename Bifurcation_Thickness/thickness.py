import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import os

# frequency = 135kHz
# Vdrc = 6.680

Vdrive = 0.1*np.arange(0,16,1)
			 #6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1
thicknesses=[[0,0,1,1,2,3,3,4,4,4,5]
			,[0,0,1,1,2,3,3,3,4,4,5]
			,[1,1,1,2,2,3,3,3,3,3,4,5]
			,[1,1,1,2,2,3,3,4,4,5,5,6,7]
			,[1,1,1,2,2,3,3,3,3,4,4,5,6,7]
			,[1,1,1,1,2,2,3,3,3,3,4,4,4,5,6]
			,[1,1,1,1,1,2,2,3,3,3,3,3,4,4,5,6]
			,[1,1,1,2,2,2,2,2,2,3,3,3,4,4,5,5,6]
			,[1,1,2,2,2,2,2,2,2,3,3,3,4,4,5,5,5,6,7]
			,[1,1,1,2,2,2,2,2,2,2,2,3,3,4,4,4,5,5,5,6,6,7]
			,[1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,4,4,4,5,5,6,6,6,7,8]
			,[1,1,2,2,2,2,2,2,2,2,2,2,2,3,3,4,4,4,4,4,4,5,5,5,6,6,6,7]
			,[1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,6,6,6,6]]
		
for i,t in zip(range(len(thicknesses)),thicknesses):
	(a,b),pcov = opt.curve_fit(lambda x,a,b: a*(x**b), 0.1*np.arange(len(t)), t)
	perr = np.sqrt(np.diag(pcov))
	print(a,perr[0],b,perr[1], abs(b-1.5247)/1.5247*100)
	if i % 6 == 0:
		plt.figure(0)
		plt.errorbar(np.arange(len(t))*0.1,t,yerr=np.ones(len(t))*0.5,marker='o',markersize=2.5,linestyle='',capsize=5)
		plt.xlabel("$V_{dr}-V_{drc}$ (V)")
		plt.ylabel("Width")
	if i == 6:
		plt.figure(1)
		plt.errorbar(np.arange(len(t))*0.1,t,yerr=np.ones(len(t))*0.5,marker='o',markersize=2.5,linestyle='',capsize=5)
		plt.plot(np.arange(len(t))*0.1,a*(np.arange(len(t))*0.1)**b)
		plt.xlabel("$V_{dr}-V_{drc}$ (V)")
		plt.ylabel("Width")
		plt.legend(["$T=3.905(V_{dr}-V_{drc})^{0.885}$","Data for $f$ = 135 kHz"])

plt.figure(0)
plt.legend(["$f$ = 105 kHz","$f$ = 135 kHz","$f$ = 165 kHz"])
plt.show()
