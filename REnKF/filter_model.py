""" data assimilation filters and state function. """

# third party imports
import numpy as np

# EnKF
def EnKF(X, HX, R, Nen):
    coeff = 1.0 / (Nen - 1.0)
    xp = X - np.mean(X, 1, keepdims=1)
    hxp = HX - np.mean(HX)
    pht = coeff * np.dot(xp, hxp)
    hpht = coeff * hxp.dot(hxp.T)
    inv = 1. / (hpht + R)
    kalman_gain_matrix = pht * inv
    return kalman_gain_matrix.reshape(2, 1)

# REnKF
def REnKF(X, HX, R, Nen):
    coeff = 1.0 / (Nen - 1.0)
    xp = X - np.mean(X, 1, keepdims=1)
    hxp = HX - np.mean(HX)
    pht = coeff * np.dot(xp, hxp)
    hpht = coeff * hxp.dot(hxp.T)
    inv = 1. / (hpht + R)
    kalman_gain_matrix = np.reshape(pht * inv, (2, 1))
    hxx = np.reshape(coeff * np.dot(hxp, xp.T), (1, 2))
    p = coeff * np.dot(xp, xp.T)
    k2_gain_matrix = kalman_gain_matrix.dot(hxx) - p
    return kalman_gain_matrix, k2_gain_matrix

# state function to propagate ensemble realizations
def state_model(theta, Nen):
    X = np.zeros((2, Nen))
    for i in range(Nen):
        X[0,i] = np.exp(-(theta[0, i] + 1)**2 - (theta[1, i] + 1)**2)
        X[1,i] = np.exp(-(theta[0, i] - 1)**2 - (theta[1, i] - 1)**2)
    return X

# observation_model
def obs_model(X):
    H = np.array([-1.5, -1.0])
    HX = H.dot(X)
    return HX
