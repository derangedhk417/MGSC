import matplotlib.pyplot as plt
import numpy             as np
import sys

if __name__ == '__main__':
	fname = sys.argv[1]

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

	plt.scatter(scales, exps, s=4)
	plt.xlabel("Scale Factor")
	plt.ylabel("Ground State Energy Prediction")
	plt.show()