import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit
from scipy.stats import norm
import math

def H(p):
	return -1*(p*np.log(p) + (1-p)*np.log(1-p))

q  = 2

blocklength = 200

generator = []
X = []
Y = []
Y1 = []
Z1 = []
Y2 = []
Z = []

generator = [1, 1]

mean  = 0.5

sigma  = 0.5

power = poly.polypow(generator, blocklength)

rv = norm(loc = blocklength*mean, scale = sigma*np.sqrt(blocklength))

q = float(q)

i = 0
for item in power:
	X.append(i)
	Y.append((item/np.power(q,  blocklength)))
	Y1.append(np.log(item)/blocklength)
	p = i/blocklength
	Z.append(H(p))
	#Y2.append(np.log(2)*(1 - (p-mean)*(p-mean)/(2*sigma*sigma*np.log(2)) - 0.5*np.log(2*3.14*blocklength*sigma*sigma)/(blocklength*np.log(2))))
	Y2.append(np.log(2)*(1 - (p-mean)*(p-mean)/(2*sigma*sigma*np.log(2))))
	#value = np.power(q, (blocklength*(1 - c/np.log(q))-0.5*np.log(2*3.14*blocklength*sigma)/np.log(q)))
	i = i+1

#plt.plot(X, Y, 'b-')
#plt.plot(X, rv.pdf(X), 'r')

plt.plot(X, Y1, 'b')
plt.plot(X, Z, 'g')
plt.plot(X, Y2, 'r')

plt.show()