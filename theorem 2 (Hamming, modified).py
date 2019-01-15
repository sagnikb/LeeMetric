'''
This program verifies theorem 2 as stated in the Hamming metric case.
The blue graph is the expected plot, and the green one is the lower bound
Clearly, the bound holds. Since we are have taken, at one point, that M = \sqrt(n), the bound does not work for small r where r is smaller than M. Can be avoided by taking r = O(n)
'''
import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit
from scipy.stats import norm
from scipy.stats import entropy
from scipy.special import binom
import math

######################################

def choose(n, k):
    return binom(n, k)
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    https://stackoverflow.com/a/3025547/10099705
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in range(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0

######################################

q = 5
alpha = int(np.floor(q/2)) #The largest value that the Lee weight can take
blocklength = 100

#####################################################
#Generates the generating Lee polynomial for given q#
#####################################################

generator = []
generator.append(1)
for i in range (1, alpha):
	generator.append(2)
if q % 2 == 0:
	generator.append(1)
else:
	generator.append(2)

#######################################
#Calculates the mean for Lee(q) metric#
#######################################

if q%2 == 0:
	mean = q/4
else:
	mean = (q*q - 1)/(4*q)

#########################################################################
#Calculates the underlying probability distribution as used in the paper#
#########################################################################

givendist = []
givendist.append(1/q)
for i in range (1, alpha):
	givendist.append(2/q)
if q % 2 == 0:
	givendist.append(1/q)
else:
	givendist.append(2/q)

###################################################
#Given p and lambda, calculates the KL divergence.# 
#Comes from the dual program for relative entropy.#
###################################################

def kllamb(p, lamb, generator):
	array = []
	number = 0
	for item in givendist:
		array.append(item*np.exp(-1*number*lamb))
		number = number + 1
	return(-1*p*lamb - np.log(np.sum(array)))

##################################################################
#Called kllamb to calculate the optimal KL divergence for given p#
##################################################################
def kl(p, generator):
	X1 = []
	Y1 = []
	for lamb in np.linspace(0, 5, 500):
		X1.append(lamb)
		Y1.append(kllamb(p, lamb, generator))
	#return(np.max(Y1))
	return(X1[np.argmax(Y1)])

############################################
#Defines f_k for the theorem 2 in the paper#
############################################

def fk(k):
    return(1/np.sqrt(choose(blocklength, k)*np.power(q-1,k)))

#########################
#q-ary entropy function #
#########################
def H(p):
    return(p*(np.log(q-1)/np.log(q)) - p*(np.log(p)/np.log(q)) - (1-p)*(np.log(1-p)/np.log(q)))

def C(nu):
    array = []
    number = 0
    for item in givendist:
        array.append(item*np.exp(-1*number*nu))
        number = number + 1
    return((np.sum(array)-array[0]))

def klnu(p, pp, nu):
    number = (p * givendist[0])/(C(nu)*(1-p))
    return(-pp*nu + p*np.log(number) - np.log(givendist[0] + number))

def kl2(p, pp):
    X = []
    Y = []
    for i in np.linspace(-5, 5, 100):
        X.append(i)
        Y.append(klnu(p, pp, i))
    if np.max(Y) > 0:
        return(np.max(Y))
    else:
        return(0)

##############
#Main Program#
##############

q = float(q)
M = int(np.sqrt(blocklength))
R = []
data = []
approx = []

for r in range(3, int(blocklength*mean - 50)):
    pp = r/blocklength
    fAf = 0
    #print(int(r/alpha), r)
    for k in range(int(r/alpha), r):
        p = k/blocklength
        fAf = fAf + np.float_power(q, blocklength - (blocklength*kl2(p,pp))/np.log(q))*(np.sqrt((q-1)*k*(blocklength-k+1)) + (q-2)*k + np.sqrt((q-1)*(k+1)*(blocklength-k)))*fk(k)*fk(k)
    ff = M
    print(r)
    R.append(r)
    data.append(fAf/ff)
    #approx.append((2*(M-1)*np.sqrt(q-1)*(np.sqrt(r*(blocklength - r)) - M)+(q-2)*(M-1)*(r - M))/M)

print(data)

plt.plot(R, data, 'b')
#plt.plot(R, approx, 'g')


plt.show()