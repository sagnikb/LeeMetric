import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit
from scipy.stats import norm
from scipy.stats import entropy
import math

###################################

blocklength = 200
q = 6

####################################

def generator(q):
	generator = []
	alpha = int(np.floor(q/2))
	generator.append(1)
	for i in range (1, alpha):
		generator.append(2)
	if q % 2 == 0:
		generator.append(1)
	else:
		generator.append(2)
	return(generator)

#####################################

def mean(q):
	if q%2 == 0:
		mean = q/4
	else:
		mean = (q*q - 1)/(4*q)
	return(mean)

######################################

def givendist(q):
	alpha = int(np.floor(q/2))
	givendist = []
	givendist.append(1/q)
	for i in range (1, alpha):
		givendist.append(2/q)
	if q % 2 == 0:
		givendist.append(1/q)
	else:
		givendist.append(2/q)
	return(givendist)

######################################

def kllamb(p, lamb, generator, q):
	array = []
	number = 0
	for item in givendist(q):
		array.append(item*np.exp(-1*number*lamb))
		number = number + 1
	return(-1*p*lamb - np.log(np.sum(array)))

######################################

def kl1(p, generator, q):
	X1 = []
	Y1 = []
	for lamb in np.linspace(0, 8, 500):
		X1.append(lamb)
		Y1.append(kllamb(p, lamb, generator, q))
	return(X1[np.argmax(Y1)])

######################################

def kl(p, generator, q):
	X1 = []
	Y1 = []
	for lamb in np.linspace(0, 8, 500):
		X1.append(lamb)
		Y1.append(kllamb(p, lamb, generator, q))
	return(np.max(Y1))

######################################
'''
X = []
Y = []
Z = []

for p in np.linspace(0.1, mean(q), 100):
	X.append(p)
	Y.append(kl1(p, generator, q))
'''
def f(p, a):
    #return(a*np.float_power(p - mean(q), 1/q))
    return(a*np.float_power(p, 1/q) - a*np.float_power(mean(q), 1/q))
'''
parameters = curve_fit(f, X, Y)

print(parameters[0][0])
'''
######################################

def kl2(p, generator, q):
    return(kllamb(p, f(p, -0.315*q - 6), generator, q))
    #return(kllamb(p, f(p, -7), generator, q))
    #return(kllamb(p, f(p, parameters[0][0]), generator, q))


######################################

X = []
Y = []
Z = []
for x in np.linspace(0, 10, 100):
    Y.append(x*x/mean(q) - 2*x + mean(q)/2)
    X.append(x)
    Z.append(0)
plt.plot(X, Y)
plt.plot(X, Z)
plt.show()