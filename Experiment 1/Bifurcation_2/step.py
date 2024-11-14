import numpy as np


for i in range(100, 20100, 100):
	uniques = set()
	time,Vdrive,VI = np.loadtxt(f'./{i}/{i}_{1:02}.csv',delimiter=',',skiprows=3).T
	uniques.update(VI)
	uniques = np.array(sorted(list(uniques)))
	if i < 3100:
		uniques = uniques/1000
	print((i<10000)*' '+f'{i/1000} V:',min(uniques[1:]-uniques[:-1]))
