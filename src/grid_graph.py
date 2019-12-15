import matplotlib.pyplot as plt
import numpy             as np
import sys

import matplotlib as mpl

if __name__ == '__main__':
	fname = sys.argv[1]


	mpl.rcParams['axes.linewidth'] = 2

	with open(fname, 'r') as file:
		raw = file.read()

	lines = raw.replace('\r\n', '\n').split('\n')

	i   = []
	j   = []
	E   = []
	err = []

	for line in lines[1:-2]:
		if not line.isspace() and line != '':
			a, b, c, d = line.replace(' ', '').split(',')
			i.append(int(a))
			j.append(int(b))
			E.append(float(c))
			err.append(float(d))

	i   = np.array(i)
	j   = np.array(j)
	E   = np.array(E)
	err = np.array(err)


	print("Minimum Energy: %f eV"%(E.min()))
	print("Max Error:      %E eV"%(err.max()))
	print("Min Error:      %E eV"%(err.min()))
	print("Mean Error:     %E eV"%(err.mean()))

	grid = np.zeros((i.max() + 1, j.max() + 1))

	for idx in range(len(i)):
		grid[i[idx]][j[idx]] = E[idx]

	N_A1 = 16
	N_A2 = 16

	A1_min = 1e5
	A1_max = 1e7
	A2_min = 1e5
	A2_max = 1e7

	fig = plt.figure(1)
	ax  = fig.add_subplot(111)
	im  = ax.imshow(grid, cmap='Blues_r')

	ax.set_aspect(1.0)

	i_space = np.linspace(0, 1, 8)
	i_label = np.sqrt(np.log(2) / np.linspace(A1_min, A1_max, 8)) * 0.2e3
	j_space = np.linspace(0, 1, 8)
	j_label = np.sqrt(np.log(2) / np.linspace(A2_min, A2_max, 8)) * 0.2e3
	ax.set_yticks(i_space * (N_A1 - 1))
	ax.set_yticklabels(["%0.3f"%i for i in i_label])
	ax.set_xticks(j_space * (N_A2 - 1))
	ax.set_xticklabels(["%0.3f"%i for i in j_label])


	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(18)

	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(18) 

	plt.xticks(rotation=70)


	cbar = fig.colorbar(im) #, shrink=0.5, aspect=20, fraction=.12,pad=.02)
	cbar.set_label(r'Energy $[\mathsf{eV}]$', size=20)
	# access to cbar tick labels:
	cbar.ax.tick_params(labelsize=20) 
	# fig.colorbar(im)
	ax.set_xlabel(r"Initial Value for Full Width Half Max (Term 2) [nm]", fontsize=20, y=-10.5)
	ax.set_ylabel(r"Initial Value for Full Width Half Max (Term 1) [nm]", fontsize=20, x=-1.5)
	plt.title(r"Nearest Local Minimum of the Expectation Value of the Hamiltonian ", fontsize=24, y=1.04)
	plt.show()