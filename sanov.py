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

def kl1(p, generator, q):
	X1 = []
	Y1 = []
	for lamb in np.linspace(0, 8, 500):
		X1.append(lamb)
		Y1.append(kllamb(p, lamb, generator, q))
	return(X1[np.argmax(Y1)])

##########################################################################################
# Same as kl1, except that it returns actual distance, not just lambda. for direct sanov #
##########################################################################################

def kl(p, generator, q):
	X1 = []
	Y1 = []
	for lamb in np.linspace(0, 8, 500):
		X1.append(lamb)
		Y1.append(kllamb(p, lamb, generator, q))
	return(np.max(Y1))

#############################################################################
# Calculation of parameters in the approximation to the Sanov theorem bound #
#############################################################################

X = []
Y = []
Z = []
for p in np.linspace(0.1, mean(q), 100):
	X.append(p)
	Y.append(kl1(p, generator, q))

# the fitting function
def f(p, a):
    return(a*np.float_power(p, 1/q) - a*np.float_power(mean(q), 1/q))

parameters = curve_fit(f, X, Y)
print(parameters[0][0])

#####################################################################################
# Uses above functions to calculate the Sanov theorem bound using the approximation #
#####################################################################################

def kl2(p, generator, q): 
    return(kllamb(p, f(p, parameters[0][0]), generator, q))

################
# Main Program # 
################

power = poly.polypow(generator(q), blocklength)

X = []
Y = []
Z = []
W = []
U = []

q = float(q) #for logs and powers to work correctly

i = 0
number = 0

for item in power[1:int(mean(q)*blocklength)]: # Plots only till peak, where the formulation holds
	if i%10 == 0: # plot only every tenth point
		X.append(i) # the x-axis
		number = number + item 
		Y.append(math.log(number, q)) # Direct Numerical Simulation
		p = i/blocklength # change to the probability interpretation
		c = (3*(4*p-q)*(4*p-q))/(2*(q*q+8)) # comes while calculating the gaussian approximation
		value = np.power(q, (blocklength*(1 - c/np.log(q))))
		value = math.log(value, q)
		Z.append(value) # gaussian approximation
		W.append(- (blocklength*(kl(p, generator, q)))/(np.log(q)) + blocklength) # direct Sanov theorem use
		U.append(- (blocklength*(kl2(p, generator, q)))/(np.log(q)) + blocklength) # sanov theorem using approximation
	i = i+1

#####################
# Plotting Commands # 
#####################

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig, ax = plt.subplots()
plt.plot(X, Y, marker = '.', color = 'b', label = r'Numerical Simulation', linewidth = '1')
plt.plot(X, Z, color = 'r', label = r'Gaussian Approximation', linewidth = '0', marker = '^', markersize = '2')
plt.plot(X, W, 'g', label = r'Full Sanov theorem approximation')
plt.plot(X, U, 'k', linestyle = 'dashed', label = r'Sanov theorem with chosen function', linewidth = '2')
plt.xlabel(r'$\displaystyle r$', fontsize = '14')
plt.ylabel(r'$\displaystyle\log_q(V_r^{(n)})$', fontsize = '14')
plt.xticks(fontsize = '12')
plt.yticks(fontsize = '12')
plt.legend(loc='lower right', title=r'$\displaystyle n = 200$, $\displaystyle q = 6$')
plt.grid(True)
fig.savefig('Sanov.svg', format='svg', dpi=1200)
plt.show()
