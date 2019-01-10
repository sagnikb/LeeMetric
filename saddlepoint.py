import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit
from scipy.stats import norm
import math

blocklength = 100

def polynomial(p,q):
	result = []
	result.append(p)
	alpha = int(np.floor(q/2))
	for i in range(1, alpha, 1):
		result.append(2*p - 2*i)
	if q%2 == 0:
		result.append(p - alpha)
	else:
		result.append(2*p - 2*alpha)
	return(result)

q = 6

alpha = int(np.floor(q/2))

generator = []
generator.append(1)
for i in range (1, alpha):
	generator.append(2)
if q % 2 == 0:
	generator.append(1)
else:
	generator.append(2)

power = poly.polypow(generator, blocklength)

X = []
Y = []
Y1 = []
Z1 = []
Z2 = []
Y2 = []
W =[]

q = float(q)

i = 1
for item in power:
	if i > 1.5*blocklength:
		break
	X.append(i)
	Y.append((item/np.power(q,  blocklength)))
	Y1.append(item)
	Y2.append(np.log(item))
	p = i/blocklength
	c = (3*(4*p-q)*(4*p-q))/(2*(q*q+8))
	value = np.power(q, (blocklength*(1 - c/np.log(q))))
	Z1.append(value)
	value = np.log(value)
	Z2.append(value)

	roots = poly.polyroots(polynomial(p, q))
	realroots = []
	for root in roots:
		if np.imag(root) == 0:
			realroots.append(np.real(root))
	if(len(realroots) == 0):
		break
	root = min(realroots)

	W.append(np.log(poly.polyval(root, power)/np.power(root, i)))

	i = i+1

#W.append(100)

plt.plot(X, Y2, 'b')
plt.plot(X, W, 'r')
plt.plot(X, Z2, 'g')
plt.show()