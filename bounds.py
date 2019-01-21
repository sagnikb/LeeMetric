import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit
from scipy.stats import norm
from scipy.stats import entropy
import math

###################################

blocklength = 200
q = 2

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

delta = []
hamming = []
gv = []
eb  = []
singleton = []

for d in np.linspace(0, 2*mean(q), 75):
    delta.append(d)
    if(d > mean(q)):
        gv.append(0)
    else:
        gv.append(kl2(d, generator, q)/np.log(q))
    hamming.append(kl2(d/2, generator, q)/np.log(q))
    if q%2==0:
        if 1-4*(q-1)*d/(q*q) > 0:
            singleton.append(1-4*(q-1)*d/(q*q))
        else:
            singleton.append(0)
    else:
        if(1-4*d/(q-1) > 0):
            singleton.append(1-4*d/(q-1))
        else:
            singleton.append(0)
    number = mean(q)*(1-np.sqrt(1-d/mean(q)))
    eb.append(kl2(number, generator, q)/np.log(q))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig, ax = plt.subplots()

#plt.plot(X, Y2, marker = '.', color = 'b', label = r'Numerical Simulation', linewidth = '1')
#plt.plot(X, Z2, color = 'r', label = r'Gaussian Approximation', linewidth = '0', marker = '^', markersize = '2')
#plt.plot(X, W2, 'g', label = r'Full Sanov theorem approximation')
#plt.plot(X, W3, 'k', linestyle = 'dashed', label = r'Sanov theorem with chosen function', linewidth = '3')

plt.plot(delta, gv, color = 'b', marker = '.', label = r'Gilbert-Varshamov bound', linewidth = '0', markersize = '3')
plt.plot(delta, hamming, 'g', marker = '^', label = r'Hamming bound', linewidth = '0', markersize = '3')
plt.plot(delta, singleton, '--k', label = r'Singleton bound', linewidth = '1', markersize = '3')
plt.plot(delta, eb, color = 'k', marker = 'v', label = r'Elias-Bassalygo bound', linewidth = '0', markersize = '3')

plt.xlabel(r'$\displaystyle \delta$', fontsize = '14')
plt.ylabel(r'$\displaystyle R(\delta)$', fontsize = '14')

plt.legend(loc='upper right', title=r'$\displaystyle q = 6$')

plt.xticks(fontsize = '12')
plt.yticks(fontsize = '12')

plt.grid(True)

#fig.savefig('myimage.pdf', format='pdf', dpi=1200)
fig.savefig('bounds.svg', format='svg', dpi=1200)
plt.show()

