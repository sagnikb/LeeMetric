import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit
from scipy.stats import norm
from scipy.stats import entropy
import math

q = 4
alpha = int(np.floor(q/2))
blocklength = 200

givendist = []
givendist.append(1/q)
for i in range (1, alpha):
	givendist.append(2/q)
if q % 2 == 0:
	givendist.append(1/q)
else:
	givendist.append(2/q)
'''
def kl(p, generator):
	#print(p)
	newdist = []
	newdist.append(1 - 2*p/(alpha + 1))
	for i in range (1, alpha):
		newdist.append(2*p/(alpha*(alpha + 1)))
	newdist.append(2*p/(alpha*(alpha + 1)))
	#print(np.sum(newdist))
	return(entropy(newdist, givendist, q))
'''

def f(p, a, b):
	return(-a*np.float_power(p, 1/b) - 4.9)

def kllamb(p, lamb, generator):
	array = []
	number = 0
	for item in givendist:
		array.append(item*np.exp(-1*number*lamb))
		number = number + 1
	#print(array)
	#print(np.log(np.sum(array)))
	return(-1*p*lamb - np.log(np.sum(array)))

def kl(p, generator):
	X1 = []
	Y1 = []
	for lamb in np.linspace(0, 5, 100):
		X1.append(lamb)
		Y1.append(kllamb(p, lamb, generator))
	#return(kllamb(p, 4*(mean-p)*(mean-p), generator))
	return(np.max(Y1))
	
	#lamb = f(p, 4.6654909, 2.95935776)
	#return(kllamb(p, -lamb, generator))

generator = []
generator.append(1)
for i in range (1, alpha):
	generator.append(2)
if q % 2 == 0:
	generator.append(1)
else:
	generator.append(2)

if q%2 == 0:
	mean = q/4
else:
	mean = (q*q - 1)/(4*q)

power = poly.polypow(generator, blocklength)

X = []
Y = []
Y1 = []
Z2 = []
Y2 = []
W2 = []

q = float(q)

i = 0
for item in power[1:int(mean*blocklength)]:
	X.append(i)
	Y.append((item/np.power(q,  blocklength)))
	Y1.append(item)
	Y2.append(math.log(item, q))
	p = i/blocklength
	c = (3*(4*p-q)*(4*p-q))/(2*(q*q+8))
	value = np.power(q, (blocklength*(1 - c/np.log(q))))
	value = math.log(value, q)
	Z2.append(value)
	W2.append(math.log(blocklength, q) - (blocklength*(kl(p, generator)))/(np.log(q)) + blocklength)
	#W2.append(math.log(blocklength, q) - (blocklength*(kl(p, generator)))/(np.log(q)))
	i = i+1
plt.plot(X, Y2, 'b')
plt.plot(X, Z2, 'r-')
plt.plot(X, W2, 'g')

plt.show()
