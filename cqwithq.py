import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit
from scipy.stats import norm
from scipy.stats import entropy
import math

###################################

blocklength = 200

#######################################################
# To generate the generating polynomial for a given q #
#######################################################
 
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

######################
# Mean for a given q #
######################

def mean(q):
	if q%2 == 0:
		mean = q/4
	else:
		mean = (q*q - 1)/(4*q)
	return(mean)

############################################
# The probability distribution for given q #
############################################

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

##########################################################################
# Calculates the Sanov's theorem D(P||Q) term using the dual formulation #
# First function calculates it for a given free param lambda             #
# Second one simply takes maximum over all lambda                        #
##########################################################################

def kllamb(p, lamb, generator, q):
	array = []
	number = 0
	for item in givendist(q):
		array.append(item*np.exp(-1*number*lamb))
		number = number + 1
	return(-1*p*lamb - np.log(np.sum(array)))

def kl(p, generator, q):
	X1 = []
	Y1 = []
	for lamb in np.linspace(0, 8, 500):
		X1.append(lamb)
		Y1.append(kllamb(p, lamb, generator, q))
	return(X1[np.argmax(Y1)])

################
# Main Program # 
################

qdata = [] #To store the q-points 
parameterdata = [] #To store c(q)

for q in range(2,200):
	X = []
	Y = []

	for p in np.linspace(0.1, mean(q), 100):
		X.append(p)
		Y.append(kl(p, generator, q))

	# The fitting function 
	def f(p, a):
		return(-a*np.float_power(p, 1/q) + a*np.float_power(mean(q), 1/q))

	parameters = curve_fit(f, X, Y)
	
	qdata.append(q)
	parameterdata.append(parameters[0][0])

	print(q, parameters[0][0])

#####################
# Plotting Commands # 
#####################

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig, ax = plt.subplots()
plt.plot(qdata, parameterdata, 'b')
plt.xlabel(r'$\displaystyle q$', fontsize = '14')
plt.ylabel(r'$\displaystyle c(q)$', fontsize = '14')
plt.xticks(fontsize = '12')
plt.yticks(fontsize = '12')
plt.grid(True)
fig.savefig('cqwithq.svg', format='svg', dpi=1200)
plt.show()
