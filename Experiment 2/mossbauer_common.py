import matplotlib.pyplot as plt
import numpy as np
import os
import scipy as sp
import scipy.optimize as opt

from lmfit.models import GaussianModel

filelist = ["alpha_Fe",		# ferro, metallic bond
			"Fe2O3",		# +3 ferro (above 260 K), coordination complex
			"K4Fe_CN_6",	# +2, Fe is coordinate
			"LiFePO4_1",	# +2, ionic
			"stainless",	# anti-ferro (due to chromium), metallic bond
			"FeC2O4_1",		# 2+, ionic
			"7mu_Fe",		# ferro, metallic bond
			"LiFePO4_2",
			"FeC2O4_2"]

latex_names = ["$\\alpha$-Fe", 
			"Fe$_2$O$_3$", 
			"K$_4$Fe(CN)$_6\\cdot$3H$_2$O", 
			"LiFePO$_4$ Set 1", 
			"Stainless Steel",
			"FeC$_2$O$_4\\cdot$2H$_2$O Set 1",
			"7$\\mu$ Fe",
			"LiFePO$_4$ Set 2",
			"FeC$_2$O$_4\\cdot$2H$_2$O Set 2"]


accepted = [np.array([-27, -16, -4.8, 3.4, 14, 25])*10,
			np.array([-37.9, -20.0, -4.5, 9.2, 25.7, 41.2])*10,
			np.array([-0.19])*10,
			np.array([-1.2, 13.1])*10,
			np.array([-1.4])*10,
			np.array([1.6, 9.9])*10,
			np.array([-27, -16, -4.8, 3.4, 14, 25])*10,
			np.array([-1.2, 13.1])*10,
			np.array([1.6, 9.9])*10]
lookup = {2:[0,0], 4:[0,1,1,0], 12:[0,1,2,3,4,5,5,4,3,2,1,0]}

exp_peaks = [12, 12, 2, 4, 2, 4, 12, 4, 4]
c = 299792458 # m/s

def file_check():
	missing_files = 0
	for f in filelist:
		if not os.path.isfile(f+".TKA"):
			print("Missing file: ", f)
		
	if missing_files: exit()
	
def v_to_E(v):
	E = 14.4e3
	return E*v/c
	
def dv(t, dt, phi, dphi):
	w = 2*np.pi
	R = 2e-3
	dR = .1e-3
	return np.sqrt((w*w*R*np.sin(w*t+phi)*dt)**2
					+(w*R*np.sin(w*t+phi)*dphi)**2
					+(w*np.cos(w*t+phi)*dR)**2)
	
def v(t, phi, model = 0):
	w = 2*np.pi
	R = 2e-3
	L = 0.1275
	if model == 0:
		return -w*R*np.cos(w*t+phi)
	return w*R*R*np.cos(w*t+phi)*np.sin(w*t+phi)/np.sqrt(L*L-R*R*np.cos(w*t+phi)*np.cos(w*t+phi))-w*R*np.cos(w*t+phi)

def sin_fit(t, c, A, phi):
	return c+A*np.sin(2*np.pi*t+phi)
	
def exact_fit(t, a, b, c, phi):
	return a - b + np.sqrt(b**2-(c*np.cos(2*np.pi*t+phi))**2)+c*np.sin(2*np.pi*t+phi)

def gauss_fit(x, a, mu, sigma):
	return a*np.exp(-(x-mu)*(x-mu)/(2*sigma*sigma))
	
def calculate_energies(result, popt, perr, index):
	result_params = list(reversed(result.best_values.keys()))
	#total_fit = np.zeros_like(bins)
	vE_data = np.zeros((exp_peaks[index], 4))
	print("t (s)\tE (neV)")
	for i in range(len(result_params)//3):
		sigma = result.best_values[result_params[int(i*3)]]
		mu = result.best_values[result_params[int(i*3+1)]]
		amp = result.best_values[result_params[int(i*3+2)]]
		vE_data[i][0] = v(mu, popt[-1])
		vE_data[i][1] = dv(mu, sigma, popt[-1], perr[-1])
		vE_data[i][2] = v_to_E(vE_data[i][0])
		vE_data[i][3] = v_to_E(vE_data[i][1])
		
		'''
		print('{:.3f}'.format(mu),
			'\t{: 4.0f}'.format(1e9*vE_data[i][2]),
			'+- {:2.0f}'.format(1e9*vE_data[i][3]),
			'\t{: 4.0f}'.format(1e9*v_to_E(v(mu, popt[-1], model=1))))
			'''
		print('{:.3f}'.format(mu),
			' & {: 4.0f}'.format(1e9*vE_data[i][2]),
			'$\\pm$ {:2.0f}'.format(1e9*vE_data[i][3]))
	print("\nPeak #\tE (neV)\t\tRef (neV)\tDiff")
	'''
	temp1 = None
	temp2 = None
	if filelist[index] == "LiFePO4_1" or filelist[index] == "LiFePO4_2":
		temp1 = np.load('LiFePO4_1_vE_data.npy').T
		temp2 = np.load('LiFePO4_2_vE_data.npy').T
		
	elif filelist[index] == "FeC2O4_1" or  filelist[index] == "FeC2O4_2":
		temp1 = np.load('FeC2O4_1_vE_data.npy').T
		temp2 = np.load('FeC2O4_2_vE_data.npy').T
	
	if temp1 is not None:
		vE_data = vE_data.T
		vE_data[2] = (temp1[2]+temp2[2])/2
		vE_data[3] = np.sqrt(temp1[3]**2+temp2[3]**2)/np.sqrt(2)
		vE_data = vE_data.T
		print(vE_data[:,2:]*1e9)
		'''
	for i in range(exp_peaks[index]//2):
		E_mu = 1e9*(vE_data[i][2]+vE_data[-1-i][2])/2
		E_sigma = 1e9*np.sqrt(vE_data[i][3]**2+vE_data[-1-i][3]**2)/np.sqrt(2)
		E_ref = accepted[index][lookup[exp_peaks[index]][i]]
		print(i+1, '\t{: 4.0f}'.format(E_mu),
			'+- {:2.0f}'.format(E_sigma),
			'\t{: 4.0f}'.format(E_ref),
			'\t{: 5.1f}'.format((E_mu-E_ref)/E_sigma))
	print()
	if exp_peaks[index] == 2:
		print("Isomer Shift:", '{:.0f}'.format(1e9*(vE_data[0][2]+vE_data[1][2])/2), '{:.0f}'.format(1e9*np.sqrt(vE_data[0][3]**2+vE_data[1][3]**2)/np.sqrt(2)))
	elif exp_peaks[index] == 4:
		E = np.zeros(2)
		dE = np.zeros(2)
		for i in range(2):
			E[i] = 1e9*(vE_data[i][2]+vE_data[-1-i][2])/2
			dE[i] = (1e9*(vE_data[i][3]+vE_data[-1-i][3])/2)**2
		print("Isomer Shift:", '{:.0f}'.format((E[0]+E[1])/2), '{:.0f}'.format(np.sqrt((dE[0]+dE[1])/2)))
		print("Quadrupole Split:", '{:.0f}'.format(E[1]-E[0]), '{:.0f}'.format(np.sqrt(dE[1]+dE[0])))
	elif exp_peaks[index] == 12:
		E = np.zeros(6)
		dE = np.zeros(6)
		for i in range(6):
			E[i] = 1e9*(vE_data[i][2]+vE_data[-1-i][2])/2
			dE[i] = (1e9*(vE_data[i][3]+vE_data[-1-i][3])/2)**2
		'''
		print("Isomer Shift:",	'{:.0f}'.format((E[0]+E[1]+E[4]+E[5])/4), '{:.0f}'.format((E[0]+E[2]+E[3]+E[5])/4),
								'{:.0f}'.format(((E[0]+E[1]+E[4]+E[5])/4 + (E[0]+E[2]+E[3]+E[5])/4)/2),
								'{:.0f}'.format(np.sqrt(((dE[0]+dE[1]+dE[4]+dE[5])/4 + (dE[0]+dE[2]+dE[3]+dE[5])/4)/2)))
		print("Quadrupole Split:",  '{:.0f}'.format((E[0]+E[5]-E[1]-E[4])/2), '{:.0f}'.format((E[0]+E[5]-E[2]-E[3])/2),
									'{:.0f}'.format(((E[0]+E[5]-E[1]-E[4])/2 + (E[0]+E[5]-E[2]-E[3])/2)/2),
									'{:.0f}'.format(np.sqrt(((dE[0]+dE[5]-dE[1]-dE[4])/2 + (dE[0]+dE[5]+dE[2]+dE[3])/2)/2)))
		print("Magnetic Split (g,1/2):", '{:.0f}'.format(E[3]-E[1]), '{:.0f}'.format(E[4]-E[2]),
										'{:.0f}'.format((E[3]-E[1]+E[4]-E[2])/2), '{:.0f}'.format(np.sqrt((dE[3]+dE[1]+dE[4]+dE[2])/2)))
		print("Magnetic Split (e,1/2):", '{:.0f}'.format(E[2]-E[1]), '{:.0f}'.format(E[4]-E[3]),
										'{:.0f}'.format((E[2]-E[1]+E[4]-E[3])/2), '{:.0f}'.format(np.sqrt((dE[2]+dE[1]+dE[4]+dE[3])/2)))
		print("Magnetic Split (e,3/2):", '{:.0f}'.format(E[5]-E[0]-(E[3]-E[1])), '{:.0f}'.format(E[5]-E[0]-(E[4]-E[2])),
										'{:.0f}'.format(E[5]-E[0]-(E[3]-E[1]+E[4]-E[2])/2), '{:.0f}'.format(np.sqrt(dE[5]+dE[0]+(dE[3]+dE[1]+dE[4]+dE[2])/2)))
		print()
		'''
		print("Isomer Shift:",  '{:.0f}'.format(((E[0]+E[1]+E[4]+E[5])/4 + (E[0]+E[2]+E[3]+E[5])/4)/2),
								'{:.0f}'.format(np.sqrt(((dE[0]+dE[1]+dE[4]+dE[5])/4 + (dE[0]+dE[2]+dE[3]+dE[5])/4)/2)))
		print("Quadrupole Split:",  '{:.0f}'.format(((E[0]+E[5]-E[1]-E[4])/2 + (E[0]+E[5]-E[2]-E[3])/2)/2),
									'{:.0f}'.format(np.sqrt(((dE[0]+dE[5]-dE[1]-dE[4])/2 + (dE[0]+dE[5]+dE[2]+dE[3])/2)/2)))
		print("Magnetic Split (g,1/2):", '{:.0f}'.format((E[3]-E[1]+E[4]-E[2])/2), '{:.0f}'.format(np.sqrt((dE[3]+dE[1]+dE[4]+dE[2])/2)))
		print("Magnetic Split (e,1/2):", '{:.0f}'.format((E[2]-E[1]+E[4]-E[3])/2), '{:.0f}'.format(np.sqrt((dE[2]+dE[1]+dE[4]+dE[3])/2)))
		print("Magnetic Split (e,3/2):", '{:.0f}'.format(E[5]-E[0]-(E[3]-E[1]+E[4]-E[2])/2), '{:.0f}'.format(np.sqrt(dE[5]+dE[0]+(dE[3]+dE[1]+dE[4]+dE[2])/2)))
		print()
	return vE_data
	
def gen_model(num_peaks, guess_centers, guess_amps):
	if num_peaks == 2:
		model = (GaussianModel(prefix='p0_')
				+GaussianModel(prefix='p1_'))

		params = model.make_params(p0_center=guess_centers[0], p0_amplitude=guess_amps[0], p0_sigma=0.002,
									p1_center=guess_centers[1], p1_amplitude=guess_amps[1], p1_sigma=0.002)
		return model, params
	elif num_peaks == 4:
		model = (GaussianModel(prefix='p0_')
				+GaussianModel(prefix='p1_')
				+GaussianModel(prefix='p2_')
				+GaussianModel(prefix='p3_'))

		params = model.make_params(p0_center=guess_centers[0], p0_amplitude=guess_amps[0], p0_sigma=0.002,
									p1_center=guess_centers[1], p1_amplitude=guess_amps[1], p1_sigma=0.002,
									p2_center=guess_centers[2], p2_amplitude=guess_amps[2], p2_sigma=0.002,
									p3_center=guess_centers[3], p3_amplitude=guess_amps[3], p3_sigma=0.002)
		return model, params
	elif num_peaks == 12:						
		model = (GaussianModel(prefix='p0_')
				+GaussianModel(prefix='p1_')
				+GaussianModel(prefix='p2_')
				+GaussianModel(prefix='p3_')
				+GaussianModel(prefix='p4_')
				+GaussianModel(prefix='p5_')
				+GaussianModel(prefix='p6_')
				+GaussianModel(prefix='p7_')
				+GaussianModel(prefix='p8_')
				+GaussianModel(prefix='p9_')
				+GaussianModel(prefix='pa_')
				+GaussianModel(prefix='pb_')) 
		params = model.make_params(p0_center=guess_centers[0], p0_amplitude=guess_amps[0], p0_sigma=0.002,
									p1_center=guess_centers[1], p1_amplitude=guess_amps[1], p1_sigma=0.002,
									p2_center=guess_centers[2], p2_amplitude=guess_amps[2], p2_sigma=0.002,
									p3_center=guess_centers[3], p3_amplitude=guess_amps[3], p3_sigma=0.002,
									p4_center=guess_centers[4], p4_amplitude=guess_amps[4], p4_sigma=0.002,
									p5_center=guess_centers[5], p5_amplitude=guess_amps[5], p5_sigma=0.002,
									p6_center=guess_centers[6], p6_amplitude=guess_amps[6], p6_sigma=0.002,
									p7_center=guess_centers[7], p7_amplitude=guess_amps[7], p7_sigma=0.002,
									p8_center=guess_centers[8], p8_amplitude=guess_amps[8], p8_sigma=0.002,
									p9_center=guess_centers[9], p9_amplitude=guess_amps[9], p9_sigma=0.002,
									pa_center=guess_centers[10], pa_amplitude=guess_amps[10], pa_sigma=0.002,
									pb_center=guess_centers[11], pb_amplitude=guess_amps[11], pb_sigma=0.002)	
		return model, params

if __name__ == "__main__":
	file_check()
	for index, f in enumerate(filelist):
		plt.subplot(3,3,index+1)
		plt.title(f)
		counts = np.loadtxt(f+".TKA")[14200:14325]
		bins = np.arange(len(counts)+1)+14200
		baseline = np.mean(counts[:50])
		quarters = [None, None]
		for i,count in enumerate(counts):
			if quarters[1] is None and count/baseline < 0.75:
				quarters[1] = i
			if quarters[0] is None and count/baseline < 0.25:
				quarters[0] = i
				break
		plt.stairs(counts/baseline, bins, color='k')
		for q in quarters:
			plt.plot([q+14200,q+14200], [0,1], color='r')
		plt.plot([14200,14325], [0.25,0.25], color='b')
		plt.plot([14200,14325], [0.75,0.75], color='b')
		print(f.ljust(12), "Uncertainty:", 70*(quarters[0]-quarters[1]), "\tus")
	plt.tight_layout
	'''
	alpha_Fe     Uncertainty: 980   us
	Fe2O3        Uncertainty: 980   us
	K4Fe_CN_6    Uncertainty: 1050  us
	LiFePO4_1    Uncertainty: 980   us
	stainless    Uncertainty: 1120  us
	FeC2O4_1     Uncertainty: 1120  us
	7mu_Fe       Uncertainty: 1120  us
	LiFePO4_2    Uncertainty: 1050  us
	FeC2O4_2     Uncertainty: 1120  us

	'''
	plt.show()