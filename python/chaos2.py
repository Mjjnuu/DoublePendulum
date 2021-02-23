import numpy as np
from math import *
import matplotlib.pyplot as plt
import doublePendulum as dp
import time

"""
The purpose of this script is to illustrate the sensitivity to initial conditions of the double 
pendulum. Both theta1 and z1 are varied by amounts specified by 'multipliers'. This results in
|multipliers|^2 different initial conditions. The (theta1,z1)-pairs of the pendulums are plotted
periodically in the same picture.
"""

dt = 0.0001
t = np.arange(0.0, 6.0, dt)

Q0 = np.array([pi/2, 0.75*pi, 0.0, 0.0])
M = [1.0, 1.0]
L = [0.2, 0.2]


multipliers = np.arange(0.99995,1.00006,0.00001)
print(multipliers)

list_of_trajectories = []

for multiplier1 in multipliers:
    for multiplier2 in multipliers:
        Q = np.array([Q0[0]*multiplier1, Q0[1], 4000*(multiplier2-1), 0.0 ])
        trajectory = dp.simulate_double_pendulum(Q,M,L,t,dt,"AB4")
        list_of_trajectories.append(trajectory)

polar_coordinate_list = []


for trajectory in list_of_trajectories:
    # Normalize the angle theta1 to interval [-pi,pi] 
    theta1 = (trajectory[0,:] + pi) % (2*pi) - pi
    z1 = trajectory[2,:]
    
    theta1_z1 = np.zeros([2,len(t)])
    
    theta1_z1[0,:] = theta1
    theta1_z1[1,:] = z1
    polar_coordinate_list.append(theta1_z1)

plt.plot(polar_coordinate_list[0][0],polar_coordinate_list[0][1],'.')
plt.show()


fps = 30.0
frame_time = 1.0 / fps
indices_in_a_frame = int(frame_time / dt)





fps = 30.0
frame_time = 1.0 / fps
indices_in_a_frame = int(frame_time / dt)

fig = plt.figure()
plt.axis([-pi, pi, -6*pi, 6*pi])
fig.show()
for time_step in range(0, len(t), indices_in_a_frame):
    # Loop over all different pendulums
    
    theta1 = []
    z1 = []
    for pendulum in range(len(list_of_trajectories)):
        theta1.append(polar_coordinate_list[pendulum][0][time_step])
        z1.append(polar_coordinate_list[pendulum][1][time_step])    
    
    fig.clear()    
    plt.plot(theta1, z1, '.')  
        
    plt.axis([-pi, pi, -6*pi, 6*pi])
    fig.canvas.draw()
    #plt.savefig('frames/theta1z1{:04d}.png'.format(time_step//indices_in_a_frame), dpi=72, bbox_inches='tight')
    plt.pause(0.0001) 
    











