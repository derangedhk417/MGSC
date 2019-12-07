from scipy.integrate import nquad
import numpy as np
import time

R   = 1.0
a_1 = 1.0
a_2 = 1.0
a_3 = 1.0

def f(r, theta, phi):
	Theta  = a_1*np.square(np.sin(theta))*np.square(np.cos(phi))
	Theta += a_2*np.square(np.sin(theta))*np.square(np.sin(phi))
	Theta += a_3*np.square(np.cos(theta))
	return np.sin(theta)*((np.square(r)*np.exp(-np.square(r)*(Theta))) * (1/np.sqrt(np.square(R) + np.square(r) - 2*r*R*np.cos(theta))))

if __name__ == '__main__':
	start = time.time_ns()
	result = nquad(f, [(0, 100), (0, np.pi), (0, 2*np.pi)], opts={"epsabs":1e-1})
	stop = time.time_ns()
	print(result)
	print("Duration: %f"%((stop - start) / 1e9))