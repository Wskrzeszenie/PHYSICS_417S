import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.optimize as opt

V_1n538 = np.array([-200,	0,	0.8,	0.9,	1.0])
I_1n538 = np.array([-0.5e-6,	0,	2,	0.7,	2.0])

popt,_ = opt.curve_fit(lambda x,a,b: a*np.exp(b*x)-0.5e-6, V_1n538, I_1n538)

L = 25e-3
L_R = 68.4
R = 220 + L_R
C = 15e-12
#def I_D(V):
#	return popt[0]*np.exp(popt[1]*V)-0.5e-6
	
def V_D(q):
	return min(0.7,q/C)
	#return np.log((I+0.5e-6)/popt[0])/popt[1]

'''	
diffeq

V_dr/2*np.cos(2*np.pi*f*t) = L*d^2q/dt^2 + R*dq/dt - V_D
'''

sz = 100000

def circuit(t, q, V_dr=0.6, f=85e3):
	return [q[1], 1/L*(-R*q[1]+V_D(q[0])+V_dr/2*np.cos(2*np.pi*f*t))]
	
t = np.linspace(0,1e-3,2*sz)

q = sp.integrate.odeint(circuit, [0, 1e-5], t, tfirst=True).T[1]
print(q)


sol_fft = sp.fft.fft(q[sz//2:])
freq = sp.fft.fftfreq(sz, 1e-3/sz)[:sz//2]

plt.figure()
plt.subplot(211)
plt.plot(t[sz//2:], R*q[sz//2:])
plt.subplot(212)
plt.plot(freq, 2.0/sz * np.abs(sol_fft[:sz//2]))
plt.xlim((0,100e3))
plt.show()
'''
sol = sp.integrate.solve_ivp(circuit, [min(t),max(t)], [-1e-6], t_eval=t)#, dense_output=True)

sol_fft = sp.fft.fft(sol.y[0].T[sz//2:])
freq = sp.fft.fftfreq(sz, 1e-3/sz)[:sz//2]

plt.figure()
plt.subplot(211)
plt.plot(t[sz//2:],sol.y[0].T[sz//2:])
plt.subplot(212)
plt.plot(freq, 2.0/sz * np.abs(sol_fft[:sz//2]))
plt.xlim((0,100e3))
plt.show()
'''
