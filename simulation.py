import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import os
import scipy as sp
import scipy.optimize as opt
import time

L = 0.025
L_R = 68.4
R = 220 + L_R
C1 = 1e-7
C2 = 1e-9
V0 = 0.05

def diode(q):
	return (C2-C1)/(2*C1*C2)*np.abs(q)+(C2+C1)/(2*C1*C2)*q+V0

def circuit(t, q, V_dr, f=15.5e3):
		return [q[1], 1/L*(-R*q[1]-diode(q[0])+V_dr/2*np.cos(2*np.pi*f*t))]

def simulate(V_dr, f=15.5e3):
	sz = 10000000*2
	t_end = 0.05
	t = np.linspace(0,t_end,sz)
	q = sp.integrate.odeint(circuit, [0, 0], t, args=(V_dr,f), tfirst=True).T
	peaks, _ = sp.signal.find_peaks(q[1]*100000, distance = (1/15.5e3)/(t_end/sz)*0.85)
	uniques = set(R*(q[1][peaks])[len(peaks)//2:])
	return uniques
	
def four_plot(V_dr, f=15.5e3, density=1):
	plt.figure(str(V_dr))
	sz = 10000000*density
	t_end = 0.05
	t = np.linspace(0,t_end,sz)
	q = sp.integrate.odeint(circuit, [0, 0], t, args=(V_dr,f), tfirst=True).T
	sol_fft = sp.fft.fft(q[1][-sz//4:])
	freq = sp.fft.fftfreq(len(q[1][-sz//4:]), t_end/sz)
	
	plt.subplot(2,2,1)
	plt.plot(t[-sz//32:],q[0][-sz//32:])
	plt.xlabel('$t$')
	plt.ylabel('$q(t)$')

	plt.subplot(2,2,2)
	plt.plot(t[-sz//32:],q[1][-sz//32:])
	plt.xlabel('$t$')
	plt.ylabel('$i(t)$')
	
	plt.subplot(2,2,3)
	plt.plot(q[0][-sz//8:],q[1][-sz//8:])
	plt.xlabel('$q(t)$')
	plt.ylabel('$i(t)$')

	plt.subplot(2,2,4)
	plt.plot(freq[:len(freq)//2], abs(sol_fft[:len(freq)//2]))
	plt.xlabel('$f$ (Hz)')
	plt.ylabel('Amplitude')
	plt.yscale('log')
	
	plt.xlim(left=0,right=f*1.1)
	plt.tight_layout()
	
def bifurcation(start, end, n):
	with mp.Pool(processes=mp.cpu_count()-1) as p:
		results = p.map(simulate, np.linspace(start, end, n))
	for vals,V in zip(results, np.linspace(start, end, n)):
		plt.scatter(np.ones(len(vals))*V, list(vals), s=0.1, color='k')
	plt.xlabel("$V_{dr}$ ($\\mathrm{V_{PP}}$)")
	plt.ylabel("$V_R$ (V)")

if __name__ == '__main__':
	plt.figure(0)
	if not os.path.isfile("full_bifurcation.npy"):
		with mp.Pool(processes=mp.cpu_count()-1) as p:
			results = p.map(simulate, np.arange(300)*0.05)
		np.save('full_bifurcation',results)
	full_bifurcation = np.load('full_bifurcation.npy',allow_pickle=True)
	for vals,V in zip(full_bifurcation, np.arange(300)*0.05):
		plt.scatter(np.ones(len(vals))*V, list(vals), s=0.1, color='k')
	plt.xlabel("$V_{dr}$ ($\\mathrm{V_{PP}}$)")
	plt.ylabel("$V_R$ (V)")
	plt.tight_layout()
	plt.savefig('full_bifurcation.png',dpi=150)
	plt.figure(1)
	if not os.path.isfile("partial_bifurcation.npy"):
		with mp.Pool(processes=mp.cpu_count()-1) as p:
			results = p.map(simulate, np.arange(400)*0.01)
		np.save('partial_bifurcation',results)
	partial_bifurcation = np.load('partial_bifurcation.npy',allow_pickle=True)
	for vals,V in zip(partial_bifurcation, np.arange(400)*0.01):
		plt.scatter(np.ones(len(vals))*V, list(vals), s=0.1, color='k')
	plt.xlabel("$V_{dr}$ ($\\mathrm{V_{PP}}$)")
	plt.ylabel("$V_R$ (V)")
	plt.tight_layout()
	plt.savefig('partial_bifurcation.png',dpi=150)
	plt.figure(2)
	if not os.path.isfile("small_bifurcation.npy"):
		with mp.Pool(processes=mp.cpu_count()-1) as p:
			results = p.map(simulate, np.arange(800)*0.0005)
		np.save('small_bifurcation',results)
	small_bifurcation = np.load('small_bifurcation.npy',allow_pickle=True)
	for vals,V in zip(small_bifurcation, np.arange(800)*0.0005):
		plt.scatter(np.ones(len(vals))*V, list(vals), s=0.1, color='k')
	plt.xlabel("$V_{dr}$ ($\\mathrm{V_{PP}}$)")
	plt.ylabel("$V_R$ (V)")
	plt.tight_layout()
	plt.savefig('small_bifurcation.png',dpi=150)
	plt.show()