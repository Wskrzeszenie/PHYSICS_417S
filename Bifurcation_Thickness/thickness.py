import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import os

# frequency = 135kHz
# Vdrc = 6.680

Vdrive = 0.1*np.arange(0,16,1)
			 #6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1
t = np.array([1,1,1,1,1,2,2,3,3,3,3,4,4,5,6,7])*0.007874
	
(a,b),pcov = opt.curve_fit(lambda x,a,b: a*(x**b), Vdrive, t)

perr = np.sqrt(np.diag(pcov))

print(perr)
print(a,b)

x = np.linspace(0, 1.5, 50)
plt.plot(x, a*(x**b))
plt.scatter(Vdrive,t)
plt.show()
