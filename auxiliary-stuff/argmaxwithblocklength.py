import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit
import math

for q in range(3,15,1):
	alpha = int(np.floor(q/2))
	generator = []
	X = []
	Y = []
	Z = []
	generator.append(1)
	for i in range (1, alpha):
		generator.append(2)
	if q % 2 == 0:
		generator.append(1)
	else:
		generator.append(2)
	for blocklength in range(0, 250, 5):
		power = poly.polypow(generator,blocklength)
		Y.append(np.argmax(power))
		X.append(blocklength)
		Z.append(np.argmax(power)/blocklength)
	#plt.plot(X, Y)
	plt.plot(X,Z)
plt.show()