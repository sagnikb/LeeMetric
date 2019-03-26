'''
This program verifies theorem 2 as stated in the Hamming metric case, and also the changed version.
The blue graph is the expected plot, and the green one is the lower bound, the red one is the one obtained by using the approximation. Clearly, red and blue coincide.
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

q = 7
alpha = int(np.floor(q/2)) #The largest value that the Lee weight can take
blocklength = 50

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

##############
#Main Program#
##############

q = float(q)
M = int(np.sqrt(blocklength))
R = []
data = []
data1 = []
approx = []
'''
for r in range(M, int(blocklength*mean)):
    fAf = 0
    for k in range(r-M+1, r):
        fAf = fAf + np.sqrt((q-1)*k*(blocklength-k+1)) + (q-2)*k + np.sqrt((q-1)*(k+1)*(blocklength-k))
    ff = M
    R.append(r)
    data.append(fAf/ff)
    approx.append((2*(M-1)*np.sqrt(q-1)*(np.sqrt(r*(blocklength - r)) - M)+(q-2)*(M-1)*(r - M))/M)
'''
'''
for r in range(M, int(blocklength*mean)):
    fAf = 0
    fAf1 = 0
    for k in range(r-M+1, r):
        fAf = fAf + np.float_power(q, blocklength*H(k/blocklength))*fk(k)*fk(k)*(np.sqrt((q-1)*k*(blocklength-k+1)) + (q-2)*k + np.sqrt((q-1)*(k+1)*(blocklength-k)))
        fAf1 = fAf1 + np.sqrt((q-1)*k*(blocklength-k+1)) + (q-2)*k + np.sqrt((q-1)*(k+1)*(blocklength-k))
    ff = M
    R.append(r)
    data.append(fAf/ff)
    data1.append(fAf1/ff)
    approx.append((2*(M-1)*np.sqrt(q-1)*(np.sqrt(r*(blocklength - r)) - M)+(q-2)*(M-1)*(r - M))/M)

plt.plot(R, data, 'r')
plt.plot(R, data1, 'b')
plt.plot(R, approx, 'g')
'''
X = []
Y = []
Z = []
W = []

for k in range(blocklength):
    X.append(k)
    Y.append(binom(blocklength, k)*np.float_power(q-1, k))
    number = np.float_power(q, blocklength*H(k/blocklength))
    Z.append(np.float_power(q, blocklength*H(k/blocklength)))
    W.append(np.float_power(q, blocklength*H(k/blocklength))/np.float_power(blocklength+1, q))
    print(k)

plt.plot(X, np.log(Y), 'r')
plt.plot(X, np.log(Z), 'b')
plt.plot(X, np.log(W), 'g')
'''
plt.plot(X, (Y), 'r')
plt.plot(X, (Z), 'b')
plt.plot(X, (W), 'g')
'''
plt.show()