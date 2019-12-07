import numpy as np
import time
from inspect import signature

def nparams(f):
	return len(signature(f).parameters)

def sign(x):
	if x >= 0:
		return 1.0
	else:
		return -1.0

def FindMinima(fn, gradfn, eta=0.1, threshold=1e-3):
	# Determine the number of arguments to the gradient.
	n_params = len(signature(gradfn).parameters)

	# Create a randomized starting position.
	position = (np.random.rand(n_params) - 0.5) * 20
	gradient = np.ones(n_params)

	while True:
		# Calculate the gradient.
		gradient = gradfn(*position)
		#print(gradient)
		
		# Check for zero gradient.
		closeness = np.abs(gradient).max()
		if closeness < threshold:
			break

		# Move each coordinate based on the gradient.
		cofactor  = ((closeness**4) / (1 + closeness**4)) + (eta / 10)
		for c in range(n_params):
			position[c] -= eta * sign(gradient[c]) * cofactor
			#position[c] -= eta * gradient[c] * cofactor
			
			
	
	# We may have reached a minimum. Return this position.
	return position

def Minimize(fn, gradfn, std_threshold=0.01, n_best=10):
	# Determine the number of arguments to the gradient.
	n_params = len(signature(gradfn).parameters)

	positions = []
	minima    = []
	# Find enough minima to start looking for convergence.
	for i in range(n_best):
		pos = FindMinima(fn, gradfn)
		val = fn(*pos)[0]
		minima.append(val)
		positions.append(pos)

	minima.sort()
	while np.array(minima[:n_best]).std() > std_threshold:
		print(len(minima))

		pos = FindMinima(fn, gradfn)
		val = fn(*pos)[0]
		minima.append(val)
		print(val)
		positions.append(pos)
		minima.sort()
	
	minidx  = np.array(minima).argmin()
	minimum = minima[minidx]
	minpos  = positions[minidx]
	return minimum, minpos

def f(x, y, z):
	return [-np.exp(-((x / 5)**2 + (y / 5)**2 + (z / 5)**2))*(np.sin(x) + np.sin(y) + np.sin(z))]

def gradf(x, y, z):
	return [
		-np.exp(-((x / 5)**2 + (y / 5)**2 + (z / 5)**2))*np.cos(x) + (2 / 25)*x*np.exp(-((x / 5)**2 + (y / 5)**2 + (z / 5)**2))*(np.sin(x) + np.sin(y) + np.sin(z)),
		-np.exp(-((x / 5)**2 + (y / 5)**2 + (z / 5)**2))*np.cos(y) + (2 / 25)*y*np.exp(-((x / 5)**2 + (y / 5)**2 + (z / 5)**2))*(np.sin(x) + np.sin(y) + np.sin(z)),
		-np.exp(-((x / 5)**2 + (y / 5)**2 + (z / 5)**2))*np.cos(z) + (2 / 25)*z*np.exp(-((x / 5)**2 + (y / 5)**2 + (z / 5)**2))*(np.sin(x) + np.sin(y) + np.sin(z))
	]



if __name__ == '__main__':
	start = time.time_ns()
	minimum, minpos = Minimize(f, gradf)
	stop = time.time_ns()
	print("Minimum: %f"%minimum)
	print("Pos:     %s"%str(minpos))
	print("Time:    %fs"%((stop - start) / 1e9))