#!/usr/bin/env python
"""
main program for parameter estimation case
"""

# standard library imports
import os
import sys
import shutil

# third party imports
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

# local imports
from utilities import read_input_data
from filter_model import EnKF, REnKF, state_model, obs_model
import penalty_fun

# convergence test for inequlity constraint
def converge_test(mean, cov, obs, sigmad, R, Nen):
    for i_test in range(50):
        X0 = np.random.multivariate_normal(mean, cov, Nen).T
        X = X0
        for iter in range(50):
            state = state_model(X, Nen)
            HX = obs_model(state)
            # misfit check
            HX_mean = np.mean(HX)
            misfit = (la.norm(obs - HX_mean))
            if (misfit > 1): break
            # get regularization parameter
            Anomaly = X - np.tile(X.mean(axis = 1), (Nen, 1)).T
            Pm = Anomaly.dot(Anomaly.T) / Nen
            lamb_m = lamb / np.linalg.norm(Pm)
            lamda = 0.5 * lamb_m * (np.tanh( 1.0 / d * (iter - S)) + 1)
            # ensemble observation
            obs_ensemble = np.random.normal(obs, sigmad, Nen).reshape(1, Nen)
            # get Kalman gain matrix
            kalman_gain_matrix, k2_gain_matrix = REnKF(X, HX, R, Nen)
            # penalty function
            penalty_list = penalty_fun.penalties(penalty_type)
            penalty_mat = np.zeros((2, Nen))
            penalty_term = np.zeros((2, Nen))
            for i in range(len(penalty_list)):
                penalty, grad_penalty = penalty_list[i](X, Nen)
                penalty_mat += lamda * W * (grad_penalty * penalty)
                penalty_term += penalty 
            # delta X
            dx1 = kalman_gain_matrix * (obs_ensemble - HX)
            dx2 = k2_gain_matrix.dot(penalty_mat)
            if iter < 10 and la.norm(dx2) > 5 * la.norm(dx1):
                dx2 = dx2 / 5
            # update X
            X = X + dx1 + dx2
        if iter < 19:
            None
        else:
            print('convergence test pass')
            break
    if i_test == 49:
        print ('convergence test fail')
    return X0

# read setups from input file
input_file = sys.argv[1]
inputs = read_input_data(input_file)
case_name = inputs['case_name'] # prior
penalty_type = int(inputs['penalty_type']) # constraint type
filter_name = inputs['filter_name']   # filter name
# hyper parameter in inflation function
lamb = float(inputs['lamb'])
S = float(inputs['S'])
d = float(inputs['d'])
Nen = int(inputs['ensemble_size'])      # ensemble size
max_iter = int(inputs['max_iteration']) # maximum itertaion number
cri = float(inputs['converge_criteria']) # convergence_criteria
plot_flag = inputs['plot_flag'] # show plots or not

# model error covariance
P = np.diag([0.01, 0.01])

# observation error variance
sigmad = 0.01
R =  sigmad * sigmad
# observation
obs = -1.0005

# weight of constraint
W = 1

# file operation
if filter_name == 'REnKF':
    folder_name = 'penalty' + str(penalty_type) + '_' + \
                filter_name + '_' + case_name + '_la' + str(lamb) \
                + '_S' + str(S) + '_d' + str(d)
else:
    folder_name = filter_name + '_' + case_name

dir = './postprocessing/' + folder_name+'/'
if os.path.exists(dir):
    shutil.rmtree('./postprocessing/' + folder_name)
    os.makedirs(dir + 'samples/theta')
    os.mkdir(dir + 'samples/HX')
else:
    os.makedirs(dir + 'samples/theta')
    os.mkdir(dir + 'samples/HX')

# prior / first guess
if case_name == 'prior2&2':
    theta_mean = np.array([2, 2])
if case_name == 'prior0&0':
    theta_mean = np.array([0, 0])
if case_name == 'prior-2&-2':
    theta_mean = np.array([-2, -2])

theta_ensemble = np.random.multivariate_normal(theta_mean, P, Nen).T

# initial variables for saving
converge_flag = 'False'
theta_all = theta_mean
HX_all = []
cost1_all = []
cost2_all = []
cost3_all = []
dx1_all = []
dx2_all = []
dx_all = []
penalty_all = []

# convergence test for inequality constraint
if filter_name == 'REnKF' and penalty_type!=0:
    theta_ensemble = converge_test(theta_mean, P, obs, sigmad, R, Nen)

for iter in range(max_iter):    
    X = theta_ensemble
    state = state_model(X, Nen)
    HX = obs_model(state)
    
    # misfit check
    HX_mean = np.mean(HX)
    HX_all.append(HX_mean)
    misfit = (la.norm(obs - HX_mean))

    # converge check
    if misfit < cri: converge_flag = 'True'
    if misfit >= 1: 
        print ('filter diverge due to statistical uncertainty, try to rerun')
        break
    # screen output
    print ('itertaion', iter, 'misfit =',misfit)
    # file operation
    with open(dir + 'results.dat', 'a') as f:
        f.write('\n'+str(iter)+' '+str(misfit)+' '+str(HX_mean) + ' ' \
                + str(np.linalg.norm(P)))
        f.close()
    np.savetxt(dir + 'samples/theta/' + 'iter' + str(iter) + '_theta.dat',\
               theta_ensemble)
    np.savetxt(dir + 'samples/HX/' + 'iter' + str(iter) + '_HX.dat', HX)

    if filter_name == 'EnKF':
        # ensemble observation
        obs_ensemble = np.random.normal(obs, sigmad, Nen).reshape(1,Nen) 
        # get Kalman gain matrix
        analysis_matrix = EnKF(X, HX, R, Nen)
        # delta X
        dx = analysis_matrix * (obs_ensemble - HX)
        # update augmented X
        X = X + dx
        # obtain the theta mean and save in theta_all
        theta_ensemble = X
        theta_mean = theta_ensemble.mean(axis = 1)
        theta_all = np.vstack((theta_all, theta_mean))
        dx = np.mean(dx, axis=1)
    else:
        # get regularization parameter
        Anomaly = X - np.tile(X.mean(axis = 1), (Nen, 1)).T
        Pm = Anomaly.dot(Anomaly.T) / Nen
        lamb_m = lamb / np.linalg.norm(Pm)
        lamda = 0.5 * lamb_m * (np.tanh( 1.0 / d * (iter - S)) + 1.0)
        # ensemble observation
        obs_ensemble = np.random.normal(obs, sigmad, Nen).reshape(1, Nen)
        # get Kalman gain matrix
        kalman_gain_matrix, k2_gain_matrix = REnKF(X, HX, R, Nen)

        # penalty function
        penalty_list = penalty_fun.penalties(penalty_type)
        penalty_mat = np.zeros((2, Nen))
        penalty_term = np.zeros((2, Nen))
        for i in range(len(penalty_list)):
            penalty, grad_penalty = penalty_list[i](theta_ensemble, Nen)
            penalty_mat += lamda * W * (grad_penalty * penalty)
            penalty_term += penalty
        # delta X
        dx1 = kalman_gain_matrix * (obs_ensemble - HX)
        dx2 = k2_gain_matrix.dot(penalty_mat)
        if iter < 10 and la.norm(dx2) > 5 * la.norm(dx1):
            print('Caution: Overcorrection')
            dx2 = dx2 / 5
        # update X
        X = X + dx1 + dx2

        # obtain theta mean and save in theta_all
        theta_ensemble = X
        theta_mean = theta_ensemble.mean(axis = 1)
        theta_all = np.vstack((theta_all, theta_mean))
        
        dx = np.mean(dx1 + dx2, axis=1)
        dx1_all.append(la.norm(dx1))
        dx2_all.append(la.norm(dx2))
        # calculate the proportion of constraint term in the cost function
        cost_3 = 0.5 * lamda * np.mean(penalty_term, axis=1).dot(
                W * np.mean(penalty_term,axis=1))
        cost3_all.append(cost_3)
        penalty_all.append(np.mean(penalty_term))
        
    # calculate cost function value of prior and model output
    coeff = 1.0 / (Nen - 1.0)
    xp = X - np.mean(X)
    p = coeff * xp.dot(xp.T)  
    cost_1 = 0.5 * np.dot(dx.T.dot(la.inv(p)), dx)
    
    HXa = obs_model(state_model(X, Nen))
    dy = np.mean(obs_ensemble - HXa)
    cost_2 = 0.5 * dy * dy / R
    
    cost1_all.append(cost_1)
    cost2_all.append(cost_2)
    
    # save the total update
    dx_all.append(la.norm(dx))
    
    if converge_flag == 'True':
        print ('reach convergence condition at iteration', iter )
        break
    if iter == (max_iter - 1): 
        print ('reach max iteration')
        iterations = np.arange(iter+1)

        if plot_flag:
            fig1 = plt.figure(1)
            plt.plot(iterations, theta_all[:max_iter, 0], label='theta1')
            plt.plot(iterations, theta_all[:max_iter, 1], label='theta2')
            plt.legend()
            plt.show()

            fig2 =plt.figure(2)
            plt.plot(iterations, cost1_all, label='cost1')
            plt.plot(iterations, cost2_all, label='cost2')

            if filter_name == 'REnKF':
                plt.plot(iterations, cost3_all, label='cost3')
                plt.legend()
                plt.show()

                fig3 =plt.figure(3)
                plt.plot(iterations, dx1_all, label='dx1')
                plt.plot(iterations, dx2_all, label='dx2')
                plt.legend()
                plt.show()
    
                fig4 = plt.figure(4)
                plt.plot(iterations, penalty_all[:max_iter], label='penalty')
                plt.legend()
                plt.show()   

# save
np.savetxt(dir + 'dx1_all.dat', dx1_all)
np.savetxt(dir + 'dx2_all.dat', dx2_all)
np.savetxt(dir + 'dx_all.dat', dx_all)
np.savetxt(dir + 'cost1_all.dat', cost1_all)
np.savetxt(dir + 'cost2_all.dat', cost2_all)
np.savetxt(dir + 'cost3_all.dat', cost3_all)
np.savetxt(dir + 'penalty_all.dat', penalty_all)
np.savetxt(dir + 'HX_all.dat', HX_all)
np.savetxt(dir + 'theta_all.dat',theta_all)
