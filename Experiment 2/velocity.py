import matplotlib.pyplot as plt
import numpy as np

def exact_v(t):
	w = 2*np.pi
	R = 2e-3
	L = 0.1275
	phi = np.pi
	return w*R*np.cos(w*t+phi)+w*R*R*np.cos(w*t+phi)*np.sin(w*t+phi)/np.sqrt(L*L-R*R*np.cos(w*t+phi)*np.cos(w*t+phi))
	
def approx_v(t):
	w = 2*np.pi
	R = 2e-3
	L = 0.1275
	phi = np.pi
	return w*R*np.cos(w*t+phi)
	
t = np.linspace(0,1,1000)
	
plt.figure(0)
plt.plot(t, (approx_v(t)-exact_v(t))*1e3)
plt.xlim(left=0,right=1)
plt.xlabel("$t (s)$")
plt.ylabel("$\\Delta_v(t)$ (mm/s)")
plt.tight_layout()
plt.savefig("velocity_plot.png", dpi=300)

plt.figure(1)
plt.plot(t, ((approx_v(t)-exact_v(t))/exact_v(t) + (approx_v(t+0.5)-exact_v(t+0.5))/exact_v(t+0.5))/2)
plt.xlim(left=0,right=0.5)
plt.xlabel("$t (s)$")
plt.ylabel("$\\frac{1}{2}\\left[\\Delta_v(t)-\\Delta_v\\left(t+\\frac{1}{2}\\right)\\right]$")
plt.tight_layout()
plt.savefig("averaged_velocity_plot.png", dpi=300)

t = 0.25
print(((approx_v(t)-exact_v(t))/exact_v(t) + (approx_v(t+0.5)-exact_v(t+0.5))/exact_v(t+0.5))/2)
print(((approx_v(t)-exact_v(t))/exact_v(t) + (approx_v(t+0.5)-exact_v(t+0.5))/exact_v(t+0.5))/2*1e-3*14.4e3/299792458*1e9)
plt.show()
