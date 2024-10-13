import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
import os

doubling = np.genfromtxt("doublings.txt",delimiter=',', filling_values=np.nan)

plt.plot(doubling[:,0],doubling[:,1:],marker='.')
plt.legend(["$V_{p2}$","$V_{p4}$","$V_{p8}$","$V_{p16}$","$V_{drc1}$","$V_{p3}$","$V_{p6}$","$V_{p12}$","$V_{drc2}$"])
plt.xlabel("Frequency $f$ (kHz)")
plt.ylabel("$V_{dr}$ ($\\mathrm{V_{PP}}$)")
plt.show()