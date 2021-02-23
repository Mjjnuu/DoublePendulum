import numpy as np
from math import *
import matplotlib.pyplot as plt
import doublePendulum as dp
import time

# Create a double pendulum with initial conditions theta1 = pi/2, theta2 = 0.75pi, z1 = 0, z2 = 0.
Q0 = np.array([pi/2, 0.75*pi, 0.0, 0.0])
# The mass of the bobs. M = [m1, m2]
M = [1.0, 1.0]
# The lenths of the pendulums: L = [l1, l2]
L = [0.2, 0.2]

# Make another double pendulum with different initial configuration
Q1 = np.array([Q0[0]*1.00001, Q0[1]*0.99999, 0.0, 0.0 ])

# The time-step
dt = 0.0001
# Simulate for 10 seconds
t = np.arange(0.0, 10.0, dt)    


# Simulate the two double pendulums and store their trajectories
trajectory_1 = dp.simulate_double_pendulum(Q0, M, L, t, dt, "AB4")
trajectory_2 = dp.simulate_double_pendulum(Q1, M, L, t, dt, "AB4")

plt.plot(t, trajectory_1[0,:])
plt.show()


# Plot the kinetic and potential energy of the first double pendulum
plt.plot(t, dp.calculate_kinetic_energy(trajectory_1, M, L))
plt.plot(t, dp.calculate_potential_energy(trajectory_1, M, L))
plt.plot(t, dp.calculate_energy(trajectory_1, M, L))
plt.show()

# Animate the pendulums
dp.animate_double_pendulum([trajectory_1, trajectory_2], L, dt)


"""
# These are for plotting the simulation times
times = [10,20,30,40,50,60]

ab4 = []
rk4 = []

for end_time in times:
    t = np.arange(0.0,end_time,dt)
    start = time.time()
    trajectory = dp.simulate_double_pendulum(Q0,M,L,t,dt,"AB4")
    end = time.time()
    
    ab4.append(end-start)
    
    start = time.time()
    trajectory = dp.simulate_double_pendulum(Q0,M,L,t,dt,"RK4")
    end = time.time()
    
    rk4.append(end-start)

print(ab4)
print(rk4)    
plt.plot(times, ab4)
plt.plot(times, rk4)
plt.show()
    
"""    
    
    
    
    
    
    
