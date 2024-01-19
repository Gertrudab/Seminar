import numpy as np 
from scipy import stats
from scipy.integrate import quad
from sklearn.neighbors import KernelDensity
import random
import math

# function to calculate log-likelihood of binned-mode
# µN - mean of N distribution
# µT - mean of T_k distribution
# t1, t2 - time interval to observed
# n_d - number of detectors

def log_likelihood (µ_n, µ_t, t1, t2, y, n_d):
    sum = 0
    for i in range(1, n_d):
        y_mean = calculate_y(i, µ_n, µ_t, t1, t2)
        sum += y[i] * math.log(y_mean) - y_mean
    return sum

# function to calculate the mean of Y_i (number of events recorded by the i_th detector)
# i - ith detector
# µN - mean of N distribution
# µT - mean of T_k distribution
# t1, t2 - time interval to observed
# r - mean number of background events
def calculate_y (i, µ_n, µ_t, t1, t2, r):
    t = t2 - t1
    return t * quad(integrand_2, t1, t2, args=(i, µ_n, µ_t, t1, t2, r))

def integrand_2(x, i, µ_n, µ_t, t1, t2, r):
    return s(i, x) * lambda_x(x=x, µ_n=µ_n, µ_t=µ_t, t1=t1, t2=t2) + r[i]

# emission density function
# x - x coordinates 
# µN - mean of N distribution
# µT - mean of T_k distribution
# t1, t2 - time interval to observed

def integrand (t, x, µ_t):
    return pdf(x, t) * 1/µ_t * math.exp(-t/µ_t)

def lambda_x (x, µ_n, µ_t, t1, t2):
    t = t2 - t1
    return µ_n/t * quad(integrand, t1 ,t2, args=(x, µ_t))


#  function that returns that joint distribution of the recorded events Y
# µN - mean of N distribution
# µT - mean of T_k distribution
# t1, t2 - time interval to observed
# n_d - number of detectors
# r - mean number of background events
def P_y (µ_n, µ_t, t1, t2, y, n_d, r):
    mult = 0
    for i in range(1, n_d):
        y_mean = calculate_y(i, µ_n, µ_t, t1, t2, r)
        mult *= 1/y[i] * math.exp(-y_mean) * y_mean ** y[i]
        
    return mult




# pdf as kernel density estimator
# x - three tuple describing the spatial location of the kth tracer atom at time t 
# function not complete: we should find t from x and give the gaussian kernel model 
# all xs recorded at the time t
def pdf (x, t):
    density = KernelDensity(kernel='gaussian', bandwidth=0.2).fit(x)
    return density.score_samples(x)

# detector unit sensitivity pattern
# s_i(x) not implemented 
def s (i, x):
    return random.uniform(0, 1)

x = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
print(pdf(x, 1))