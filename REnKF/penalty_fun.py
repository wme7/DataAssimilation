"""
constraint function
"""

import numpy as np

# theta1 + theta2 - 2 = 0
def penalty_0(theta_ensemble, Nen):
    penalty_0 = theta_ensemble[0, :] + theta_ensemble[1, :] - 2
    grad_penalty_0 = np.ones((2, Nen))
    return penalty_0.reshape(1, Nen), grad_penalty_0


# - theta1 - theta2 + 1 < 0
def penalty_1(theta_ensemble, Nen):
    penalty_1 = np.zeros((1, Nen))
    grad_penalty_1 = np.zeros((2, Nen))
    constraint = - theta_ensemble[0, :] - theta_ensemble[1, :] + 1
    for i in range(Nen):
        if constraint[i] >= 0.0:
            penalty_1[0, i] = constraint[i]**2
            grad_penalty_1[:, i] = -2 * constraint[i] * np.ones(2)
        else:
            penalty_1[0, i] = 0
            grad_penalty_1[:, i] = np.array([0, 0])
    return penalty_1, grad_penalty_1

# theta1 + theta2 - 3 < 0
def penalty_12(theta_ensemble, Nen):
    penalty_12 = np.zeros((1, Nen))
    grad_penalty_12 = np.zeros((2, Nen))
    constraint = theta_ensemble[0, :] + theta_ensemble[1, :] - 3
    for i in range(Nen):
        if constraint[i] >= 0.0:
            penalty_12[0, i] = constraint[i]**2
            grad_penalty_12[:, i] = 2 * constraint[i] * np.ones(2)
        else:
            penalty_12[0, i] = 0
            grad_penalty_12[:, i] = np.array([0, 0])
    return penalty_12, grad_penalty_12

def penalties(ind):
    if ind == 0: penalty_list = [penalty_0]
    if ind == 1: penalty_list = [penalty_1]
    if ind == 2: penalty_list = [penalty_1, penalty_12]
    return penalty_list  
