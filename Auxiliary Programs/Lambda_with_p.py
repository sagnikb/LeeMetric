import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit
from scipy.stats import norm
from scipy.stats import entropy
import math

###################################

q = 5
alpha = int(np.floor(q/2))
blocklength = 200

####################################

generator = []
generator.append(1)
for i in range (1, alpha):
	generator.append(2)
if q % 2 == 0:
	generator.append(1)
else:
	generator.append(2)

#####################################

if q%2 == 0:
	mean = q/4
else:
	mean = (q*q - 1)/(4*q)

######################################

givendist = []
givendist.append(1/q)
for i in range (1, alpha):
	givendist.append(2/q)
if q % 2 == 0:
	givendist.append(1/q)
else:
	givendist.append(2/q)

######################################

def kllamb(p, lamb, generator):
	array = []
	number = 0
	for item in givendist:
		array.append(item*np.exp(-1*number*lamb))
		number = number + 1
	return(-1*p*lamb - np.log(np.sum(array)))

def kl(p, generator):
	X1 = []
	Y1 = []
	for lamb in np.linspace(0, 8, 1000):
		X1.append(lamb)
		Y1.append(kllamb(p, lamb, generator))
	#return(np.max(Y1))
	return(X1[np.argmax(Y1)])

def f(p, a):
	return(a*np.float_power(p, 1/q) - a*np.float_power(mean, 1/q))
	return(a*np.float_power(p, 1/b) - a*np.float_power(mean, 1/b))

X = []
Y = []
Z = []

for p in np.linspace(0.1, mean, 100):
	X.append(p)
	Y.append(kl(p, generator))
	#Z.append(4*(mean-p)*(mean-p))
	#Z.append(f(p, parameters[0][0], parameters[0][1])+0.1)
	#Y.append(np.log(kl(p,generator)))
	print(p)
'''
def f(p, a):
	return(a*(p-mean)*(p-mean))
	#return(b*np.exp(-p*a) - b*np.exp(-mean*a))
'''
parameters = curve_fit(f, X, Y)

print(parameters[0][0])

for p in X:
	Z.append(f(p, parameters[0][0]))

plt.plot(X, Y, 'b')
plt.plot(X, Z, 'g')
plt.show()