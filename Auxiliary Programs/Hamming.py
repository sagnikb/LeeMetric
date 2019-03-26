import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit
import math

maxblocklength = 400

###########################
### THE HAMMING METRIC ####
###########################


###Storage for full data
data = []

generator = [1, 1]
pol = generator

for blocklength in range(1, maxblocklength):
	bldata = []
	pol = poly.polymul(pol, generator)
	for p in np.linspace(0, 0.5,100):
		if math.floor(p*(blocklength+1)) != 0:
			pdata = []
			pdata.append(blocklength + 1)
			pdata.append(p)
			pdata.append(math.floor(p*(blocklength+1)))
			pdata.append(np.max(pol[:math.floor(p*(blocklength+1))]))
			pdata.append(sum(pol[:math.floor(p*(blocklength+1))]))
			pdata.append(np.argmax(pol[:math.floor(p*(blocklength+1))]))
			bldata.append(pdata)
			#print(pdata)
	data.append(bldata)

X = []
Y = []
Z = []

for item in data[200]:
	X.append(item[1])
	Y.append((np.log(item[3]))/(item[0]))
	Z.append((np.log(item[4]))/(item[0]))

W = []

for value in X:
	w = value*np.log(value) + (1-value)*np.log(1-value)
	W.append(-1*w)

plt.plot(X,Y, 'y-')
plt.plot(X,Z, 'r-')
plt.plot(X,W, 'b-')
plt.show()