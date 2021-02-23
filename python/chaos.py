import numpy as np
from math import *
import matplotlib.pyplot as plt
import doublePendulum as dp
import time

"""
The purpose of this script is to illustrate the sensitivity to initial conditions of the double 
pendulum. Both theta1 and theta2 are varied by amounts specified by 'multipliers'. This results in
|multipliers|^2 different initial conditions. The (theta1,theta2)-pairs of the pendulums are plotted
periodically in the same picture.
"""

dt = 0.0001
t = np.arange(0.0, 6.0, dt)

# Initial condition of a pendulum. This is later varied
Q0 = np.array([pi/2, 0.75*pi, 0.0, 0.0])
M = [1.0, 1.0]
L = [0.2, 0.2]


multipliers = np.arange(0.99999,1.00002,0.00001)
print(multipliers)

# Store the trajectories of all the pendulums.
list_of_trajectories = []

for multiplier1 in multipliers:
    for multiplier2 in multipliers:
        # Change both theta1 and theta2 of each pendulum
        Q = np.array([Q0[0]*multiplier1, Q0[1]*multiplier2, 0.0, 0.0 ])
        # Simulate each double pendulum and store the trajectory in list_of_trajectories
        trajectory = dp.simulate_double_pendulum(Q,M,L,t,dt,"AB4")
        list_of_trajectories.append(trajectory)

polar_coordinate_list = []

for trajectory in list_of_trajectories:
    # Normalize the angles theta1 and theta2 to interval [-pi,pi] by mapping interval [0,pi) onto 
    # itself using a identity mapping, and the interval [pi,2pi) onto [-pi,0).
    theta1_theta2 = (trajectory[0:2,:] + pi) % (2*pi) - pi
    polar_coordinate_list.append(theta1_theta2)

# Plot the (theta1,theta2) time-line of the first double pendulum 
plt.plot(polar_coordinate_list[0][0],polar_coordinate_list[0][1],'.')
plt.show()

"""
euclidean_trajectory_list = []
for trajectory in list_of_trajectories:
    euclidean_trajectory_list.append(dp.polar_to_euclidean(trajectory,L))
"""   

fps = 30.0
frame_time = 1.0 / fps
indices_in_a_frame = int(frame_time / dt)

draw_radius = L[0] + L[1]
    
# Create a new figure
fig = plt.figure()
# Set the axes
plt.axis([-draw_radius, draw_radius, -draw_radius, draw_radius])
plt.gca().set_aspect('equal', adjustable='box')
plt.axis('off')
fig.show()
    
"""    
# These are for plotting pendulums
# For every indices_in_a_frame:th index draw the configuration of the pendulum
for i in range(0, len(list_of_trajectories[0][0]), indices_in_a_frame):
    start = time.time()
    fig.clear()
    for trajectory in euclidean_trajectory_list:
        [X1,Y1,X2,Y2] = trajectory
        plt.plot([0.0, X1[i], X2[i]], [ 0.0, Y1[i], Y2[i] ], c='#016EEE', ls='-', marker='o')
    
    
    plt.axis([-draw_radius, draw_radius, -draw_radius, draw_radius])
    plt.gca().set_aspect('equal', adjustable='box')
    plt.axis('off')
    fig.canvas.draw()
    #plt.savefig('frames/_img{:04d}.png'.format(i//frame_index), dpi=72)
    #plt.savefig('frames/image{:d}.png'.format(i//indices_in_a_frame), dpi=72, bbox_inches='tight')
    #plt.savefig('frames/chaos{:04d}.png'.format(i//indices_in_a_frame), dpi=72, bbox_inches='tight')
    end = time.time()
    
    plt.pause(0.00001)
"""


# These are for plotting the changes in (theta1,theta2)-pairs
fps = 30.0
frame_time = 1.0 / fps
indices_in_a_frame = int(frame_time / dt)

fig = plt.figure()
# Plot on the 2pi x 2pi square centered at the origin
plt.axis([-pi, pi, -pi, pi]) 
plt.gca().set_aspect('equal', adjustable='box')
fig.show()
# Loop over all time-steps in t
for time_step in range(0, len(t), indices_in_a_frame):
    # Loop over all different pendulums and store the theta1 and theta2 values of each in these 
    # vectors
    theta1 = []
    theta2 = []
    for pendulum in range(len(list_of_trajectories)):
        theta1.append(polar_coordinate_list[pendulum][0][time_step])
        theta2.append(polar_coordinate_list[pendulum][1][time_step])
    fig.clear()    
    plt.plot(theta1, theta2, '.')  
    plt.axis([-pi, pi, -pi, pi])
    plt.gca().set_aspect('equal', adjustable='box') 
    fig.canvas.draw()
    #plt.savefig('frames/thetas{:04d}.png'.format(time_step//indices_in_a_frame), dpi=72, bbox_inches='tight')
    plt.pause(0.0001) 

















