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
	for lamb in np.linspace(0, 5, 500):
		X1.append(lamb)
		Y1.append(kllamb(p, lamb, generator))
	#return(np.max(Y1))
	return(X1[np.argmax(Y1)])

X = []
Y = []
Z = []
'''
def f(p, a, b):
	return(a*np.float_power(p, 1/b) - 5)

for p in np.linspace(0.1, mean, 1000):
	X.append(p)
	Y.append(-kl(p, generator))
	#Y.append(np.log(kl(p,generator)))

parameters = curve_fit(f, X, Y)
'''
'''
for p in X:
	Z.append(f(p, parameters[0][0], parameters[0][1]))
'''
X = []
Y = []
Z = []

for p in np.linspace(0, mean, 1000):
	X.append(p)
	Y.append(kl(p, generator))
	#Z.append(4*(mean-p)*(mean-p))
	#Z.append(f(p, parameters[0][0], parameters[0][1])+0.1)
	#Y.append(np.log(kl(p,generator)))


#print(parameters[0])

plt.plot(X, Y, 'b')
#plt.plot(X, Z, 'g')
plt.show()