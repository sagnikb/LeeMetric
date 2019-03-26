import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit
from scipy.stats import norm
import math

q  = 6
alpha = int(np.floor(q/2))
blocklength = 200

if q%2 == 0:
	mean = q/4
else:
	mean = (q*q - 1)/(4*q)

if q%2 == 0:
	sigma = 0.25*np.sqrt((q*q + 8)/3)
else:
	sigma = 0.25*(1/q)*np.sqrt(((q-1)*(q+1)*(q*q + 3))/3)

generator = []
X = []
Y = []
Y1 = []
Z1 = []
Z2 = []
Y2 = []
generator.append(1)
for i in range (1, alpha):
	generator.append(2)
if q % 2 == 0:
	generator.append(1)
else:
	generator.append(2)

power = poly.polypow(generator, blocklength)

rv = norm(loc = blocklength*mean, scale = sigma*np.sqrt(blocklength))

q = float(q)

i = 0
for item in power[1:int(mean*blocklength)]:
	X.append(i)
	Y.append((item/np.power(q,  blocklength)))
	Y1.append(item)
	Y2.append(np.log(item))
	p = i/blocklength
	c = (3*(4*p-q)*(4*p-q))/(2*(q*q+8))
	value = np.power(q, (blocklength*(1 - c/np.log(q))))
	#value = np.power(q, (blocklength*(1 - c/np.log(q))-0.5*np.log(2*3.14*blocklength*sigma)/np.log(q)))
	#value = value*np.sqrt(24/(3.14*blocklength*(q*q + 8)))
	Z1.append(value)
	value = np.log(value)
	Z2.append(value)
	i = i+1

#print(Y)

#plt.plot(X, Y1, 'b-')
#plt.plot(X, Z1, 'r-')

plt.plot(X, Y2, 'b-')
plt.plot(X, Z2, 'r-')
#plt.plot(X, rv.pdf(X), 'r')

plt.show()