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

	scales = []
	errors = []
	exps   = []

	for line in lines[1:]:
		if not line.isspace() and line != '':
			a, b, c = line.replace(' ', '').split(',')
			scales.append(float(a))
			errors.append(float(b))
			exps.append(float(c))

	scales = np.array(scales)
	errors = np.array(errors)
	exps   = np.array(exps)

	print("Minimum Energy: %f eV"%(exps.min()))
	print("a:              %E"%(scales[exps.argmin()]))
	print("Max Error:      %E eV"%(errors.max()))
	print("Min Error:      %E eV"%(errors.min()))
	print("Mean Error:     %E eV"%(errors.mean()))

	scales = np.sqrt(np.log(2) / scales) * 0.2e3
	mask   = scales < 0.428
	scales = scales[mask]
	exps   = exps[mask]

	fig = plt.figure(1)
	ax  = fig.add_subplot(111)

	buffer_x = (max(scales) - min(scales)) * 0.05
	buffer_y = (max(exps)   - min(exps))   * 0.05

	_domain = (
		min(scales) - buffer_x,
		max(scales) + buffer_x
	)

	_range = (
		min(exps) - buffer_y,
		max(exps) + buffer_y
	)

	major_ticks_x = np.linspace(
		_domain[0], _domain[1],
		10
	)

	major_ticks_y = np.linspace(
		_range[0], _range[1],
		10
	)

	ax.set_xticks(major_ticks_x)
	ax.set_yticks(major_ticks_y)
	ax.grid(which='both', linestyle=':', linewidth=2)
	ax.grid(which='major', alpha=0.9)

	plt.xlim(_domain)
	plt.ylim(_range)

	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(18)

	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(18) 


	
	ax.scatter(scales[::100], exps[::100], s=16)
	ax.set_xlabel(r"Gaussian Full Width Half Max $[\mathsf{nm}]$", fontsize=20)
	ax.set_ylabel(r"$E_{gs}\;\mathsf{[eV]}$", fontsize=20)
	plt.title("Expectation Value of the Hamiltonian as a function of Gaussian Full Width Half Max", fontsize=24)
	plt.show()