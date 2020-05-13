#!/usr/bin/env python
"""
plot the contour of the posterior
"""

# standard library imports
import os
import sys

# third party imports
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpa
import matplotlib.lines as mlines

# local imports
from utilities import read_input_data

dir = './figures'
if os.path.exists(dir):
    None
else:
    os.mkdir(dir)

case_list = ['prior2&2', 'prior0&0', 'prior-2&-2']

input_file = sys.argv[1]
inputs = read_input_data(input_file)
filter_name = inputs['filter_name']   # filter name
penalty_type = int(inputs['penalty_type']) # constraint type

# hyper parameter in inflation function
lamb = float(inputs['lamb'])
S = float(inputs['S'])
d = float(inputs['d'])

params = {
        'axes.labelsize': 20,
        'font.size': 20,
        'legend.fontsize': 12,
        'xtick.labelsize': 20,
        'ytick.labelsize': 20,
        'font.family': 'serif',
        } 

matplotlib.rcParams.update(params)
fig, ax = plt.subplots(figsize=(8,7))
fig.subplots_adjust(bottom=0.2) 
ax.set_xlabel(r'$\omega_1$')
ax.set_ylabel(r'$\omega_2$')
plt.axis('equal')

ax.set_xlim([-3, 3])
ax.set_ylim([-3, 3])

# contour plot
theta1 = np.arange(-3, 4, 0.01)
theta2 = np.arange(-3, 4, 0.01)
N = len(theta1)
X, Y = np.meshgrid(theta1, theta2)

for i in range(N):
    Z1 = np.exp(-(X + 1)**2 - (Y + 1)**2)
    Z2 = np.exp(-(X - 1)**2 - (Y - 1)**2)

HX = -1.5 * Z1 - 1.0 * Z2
y = -1.0005
I_theta = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        I_theta[i,j] = np.linalg.norm(HX[i,j] - y) **2
levels = np.array([0.0, 0.1, 0.2, 0.4, 0.6, 0.8])
contour = plt.contour(X, Y, I_theta, levels)
ax.clabel(contour, inline=1, fmt='%1.1f', fontsize=16)

group1 = np.ones(2)
plt.plot(group1[0], group1[1], 'r+', markersize=12)

circle = mpa.Circle((-1,-1),np.sqrt(np.log(1.5)), edgecolor = 'r', 
                    linestyle='-', lw=3.0, fill=False)
ax.add_patch(circle)

for case_name in case_list:
    if filter_name == 'REnKF':
        folder_name = 'penalty' + str(penalty_type) + '_' + filter_name + '_' \
                        + case_name + '_la' + str(lamb)+ '_S' + str(S) + '_d' \
                        + str(d)
    else:
        folder_name = filter_name + '_' + case_name

    dir = './postprocessing/' + folder_name

    theta_all = np.loadtxt(dir + '/theta_all.dat')

    theta_prior = theta_all[0, :]
    theta_posterior = theta_all[-1, :]
    print(str(theta_posterior))
    plt.plot(theta_prior[0], theta_prior[1], 'go', markersize=12)
    plt.plot(theta_posterior[0], theta_posterior[1], 'b>', markersize=12)

        
    style="Simple,tail_width=0.5,head_width=4,head_length=8"

    kw = dict(arrowstyle=style, color="k")
    if case_name == 'prior-2&-2':
        arrow = mpa.FancyArrowPatch(theta_prior, theta_posterior, 
                             connectionstyle="arc3,rad=-0.5", **kw)
    else:
        arrow = mpa.FancyArrowPatch(theta_prior, theta_posterior, 
                             connectionstyle="arc3,rad=-.5", **kw)
    ax.add_patch(arrow)


if filter_name=='EnKF':
    penalty_type == 0
elif penalty_type == 0:
    plt.plot(theta1, - theta1 + 2, 'k-')
elif penalty_type == 1:
    y1 = theta1
    y2 = - theta1+1
    y3 = -theta1+ 10
    plt.plot(y1, y2, 'k-')
    ax.fill_between(y1, y2, y3, facecolor='blue',alpha=0.3)
elif penalty_type == 2:
    y1 = theta1
    y2 = - theta1+1
    y3 = - theta1+3
    plt.plot(y1, y2, 'k-')
    plt.plot(y1, y3, 'k-')
    ax.fill_between(y1,y2,y3, facecolor='blue',alpha=0.3)

 
# plot legend
prior = mlines.Line2D([], [], color='green', marker='o', linestyle='None',
                          markersize=12, label='prior mean')
posterior = mlines.Line2D([], [], color='blue', marker='>', linestyle='None',
                          markersize=12, label='posterior mean')

local = mlines.Line2D([], [], color='red', linestyle='-',
                          lw=3, label='local minimums')

truth = mlines.Line2D([], [], color='red', marker='+', linestyle='None',
                          markersize=12, label='truth')


plt.legend(handles=[truth, local, prior, posterior], loc='upper center', 
          bbox_to_anchor=(0.48, 1.15), ncol=4)
#plt.tight_layout()

plt.savefig('./figures/'+ filter_name+ '_penalty'+ str(penalty_type) + '_la'
            + str(lamb) + '_S' + str(S) + '_d' + str(d) + '_contour.pdf')
