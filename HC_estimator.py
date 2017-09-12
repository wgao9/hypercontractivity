#Copyright Weihao Gao, UIUC

from math import log,pi,exp,sqrt
import numpy.random as nr
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

#Main Usage Function
def HC(x, y, bandwidth=1.06, n_trial = 10, n_iter=100, sigma = 0.1, eta = 0.1):

	'''
		Estimating Hypercontractivity s(X;Y) from samples (x_i, y_i)
		See section 3 of arxiv.org/abs/1709.03XXX for details

		Input: 
		--x: 2D list of size N*d_x (or 1D list of size N if d_x = 1)
		--y: 2D list of size N*d_y (or 1D list of size N if d_y = 1)
		
		Output: 
		--One scalar s(X;Y)

		Parameters:
		--bandwidth: constant in the rule-of-thumb bandwidth selection, 
			     i.e., bw = bandwidth*std(data)*n**(-1/(d+4))
		--n_trial: number of initializaations to try to achieve the global maximum
		--n_iter: number of iteration in gradient descent
		--delta: Initialization of gradient is 1+N(0,sigma**2)
		--eta: step size in gradient descent

	'''

	assert len(x)==len(y), "Lists should have same length"
	n = len(x)
	if x.ndim == 1:
		x = x.reshape((n,1))
   	dx = len(x[0])   	
	if y.ndim == 1:
		y = y.reshape((n,1))
	dy = len(y[0])

	# Compute the bandwidth for KDE using rule-of-thumb bandwidth selection
	bw_x = bandwidth*np.std(x)*n**(-1.0/(dx+4)) 
	bw_y = bandwidth*np.std(y)*n**(-1.0/(dy+4)) 

	# Get the matrix A(j,i) = P_{XY}(x_i,y_j)/P_X(x_i)P_Y(y_j) using KDE 
	# and normalize it such that A is doubly stochastic
	A = get_PMI(x,y,bw_x,bw_y)
	A = doubly_stochastic_normalize(A)
	
	# Using gradient descent to solve the optimization problem
	# Try n_trial different initialization and take the max
	s_xy = np.zeros(n_trial)
	for T in range(n_trial):
		weight = (np.ones(n) + sigma*nr.normal(0,1,n)).clip(1e-8,np.sqrt(n))
        	weight = weight/np.mean(weight)

		for i in range(n_iter):
			obj, grad = compute(A, weight)
			weight += eta*np.sqrt(n)*grad
			weight = weight.clip(1e-8,np.sqrt(n))
			weight = weight/np.mean(weight)
		s_xy[T] = exp(obj)
	return max(s_xy)


#Compute the matrix A(j,i) = P_{XY}(x_i,y_j)/P_X(x_i)P_Y(y_j) using KDE
def get_PMI(x,y,bw_x,bw_y):

	n = len(x)
	Wx, Wy = np.zeros((n,n)), np.zeros((n,n))
	for i in range(n):
		for j in range(n):
			Wx[i][j] = exp(-la.norm(x[i]-x[j])**2/(2*bw_x**2))
			Wy[i][j] = exp(-la.norm(y[i]-y[j])**2/(2*bw_y**2))
		Wx[i] = Wx[i]/sum(Wx[i])
		Wy[i] = Wy[i]/sum(Wy[i])
	A = np.dot(Wy, Wx.transpose())
	return A

# Doubly Stochastic normalization of a matrix 
def doubly_stochastic_normalize(X):
	n = len(X)
	return X + (1.0/n+sum(sum(X))/n**2)*np.ones((n,n)) - (np.dot(X,np.ones((n,n))) + np.dot(np.ones((n,n)),X))/n

def f(x):
	return x*log(x)

def g(x):
	return 1+log(x)

#Evaluate objective function and gradient given A and w
def compute(A, w):
	n = len(w)
	v = np.dot(A,w).clip(1e-8,np.sqrt(n))
	Dx, Dy = sum(map(f,w)), sum(map(f,v))
	obj = log(Dy) - log(Dx)
	grad = np.dot(A.transpose(),map(g,v))
	grad = grad/Dy - map(g,w)/Dx

	return obj, (grad-np.mean(grad))/la.norm(grad-np.mean(grad))



