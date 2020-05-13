#!/usr/bin/env python
"""
post-process the analyzed data:
    plot evolution of parameter, dx and cost function value with respect 
    to iterations
"""

# standard library imports
import os
import sys

# third party imports
import numpy as np
import matplotlib.pyplot as plt

# local imports
from utilities import read_input_data

# to determine the iteration number
N = 50

# read setups from input file
input_file = sys.argv[1]
inputs = read_input_data(input_file)
penalty_type = int(inputs['penalty_type']) # constraint type
filter_name = inputs['filter_name']   # filter name
# hyper parameter in inflation function
lamb = float(inputs['lamb'])
S = float(inputs['S'])
d = float(inputs['d'])

fig1 = plt.figure(1, figsize=(12,6))
dir = './figures'
if os.path.exists(dir):
    None
else:
    os.mkdir(dir)

for i in range(3):
    if i == 0: case_name = 'prior-2&-2'
    if i == 1: case_name = 'prior0&0'
    if i == 2: case_name = 'prior2&2'

    if filter_name == 'REnKF':
        folder_name = 'penalty' + str(penalty_type) + '_' + filter_name + '_' \
        + case_name + '_la' + str(lamb) + '_S' + str(S) + '_d' + str(d)
    else:
        folder_name = filter_name + '_' + case_name

    dir = './postprocessing/' + folder_name

    theta_all = np.loadtxt(dir + '/theta_all.dat')
    if filter_name == 'REnKF':
        dx1_all = np.loadtxt(dir + '/dx1_all.dat')
        dx2_all = np.loadtxt(dir + '/dx2_all.dat')
        cost1_all = np.loadtxt(dir + '/cost1_all.dat')
        cost2_all = np.loadtxt(dir + '/cost2_all.dat')
        cost3_all = np.loadtxt(dir + '/cost3_all.dat')
    # theta vs iterations
    iterations = np.arange(N)
    ax1 = fig1.add_subplot(1,3,i+1)    
    ax1.set_xlabel('iterations',fontsize=15)
    ax1.set_ylabel(r'$\omega$', fontsize=15)
    ax1.set_title(case_name, fontsize=15)
    plt.plot(iterations, theta_all[:N, 0], 'b-', lw = 1.5, label=r'$\omega_1$')
    plt.plot(iterations, theta_all[:N, 1], 'r-', lw = 1.5, label=r'$\omega_2$')
    plt.legend()
    plt.tight_layout()
plt.savefig('./figures/' + 'penalty' + str(penalty_type) + '_' + filter_name +\
            '_' + 'converge.pdf')

if filter_name == 'REnKF':
    fig2 = plt.figure(2, figsize=(12,6))
    for i in range(3):
        if i == 0: case_name = 'prior-2&-2'
        if i == 1: case_name = 'prior0&0'
        if i == 2: case_name = 'prior2&2'

        dir = './postprocessing/' + folder_name

        dx1_all = np.loadtxt(dir + '/dx1_all.dat')
        dx2_all = np.loadtxt(dir + '/dx2_all.dat')
        # theta vs iterations
        iterations = np.arange(N)
        ax2 = fig2.add_subplot(1,3,i+1)
        ax2.set_xlabel('iterations',fontsize=15)
        ax2.set_ylabel(r'$dx$', fontsize=15)
        ax2.set_title(case_name, fontsize=15)

        plt.plot(iterations, dx1_all[:N], 'b-', lw = 1.5, label=r'$dx_1$')
        plt.plot(iterations, dx2_all[:N], 'r-', lw = 1.5, label=r'$dx_2$')
        plt.legend()
        plt.tight_layout()

    plt.savefig('./figures/' + 'penalty' + str(penalty_type) + '_' + 
                filter_name + '_' + 'dx.pdf')
    
    fig3 = plt.figure(3, figsize=(12,6))
    for i in range(3):
        if i == 0: case_name = 'prior-2&-2'
        if i == 1: case_name = 'prior0&0'
        if i == 2: case_name = 'prior2&2'

        dir = './postprocessing/' + folder_name

        cost1_all = np.loadtxt(dir + '/cost1_all.dat')
        cost2_all = np.loadtxt(dir + '/cost2_all.dat')
        cost3_all = np.loadtxt(dir + '/cost3_all.dat')
        # theta vs iterations
        iterations = np.arange(N)
        ax3 = fig3.add_subplot(1,3,i+1)
        ax3.set_xlabel('iterations',fontsize=15)
        ax3.set_ylabel('cost function', fontsize=15)
        ax3.set_title(case_name, fontsize=15)
        plt.plot(iterations, cost1_all[:N], 'b-', lw = 1.5, label=r'cost 1')
        plt.plot(iterations, cost2_all[:N], 'r-', lw = 1.5, label=r'cost 2')
        plt.plot(iterations, cost3_all[:N], 'g-', lw = 1.5, label=r'cost 3')
        plt.legend()
        plt.tight_layout()
    plt.savefig('./figures/' + 'penalty' + str(penalty_type) + '_' + 
                filter_name + '_' + 'cost.pdf')


