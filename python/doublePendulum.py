import matplotlib.pyplot as plt
import numpy as np
from math import *
import time

g = 9.81

"""
Calculates the Euclidean coordinates of the bobs of the double pendulum.

Input:
Q: The state vector of the system. Q[0] = theta1, Q[1] = theta2, Q[2] = z1 and Q[3] = z2
L: A vector containing the lengths of the pedulums; L = [l1, l2]

Output:
xn: The x-coordinate of the n:th bob.
yn: The y-coordinate of the n:th bob.
"""
def polar_to_euclidean(Q, L):
    [l1, l2] = L
    [theta1, theta2] = Q[0:2,:]
    # Apply the appropriate transformations. Use np.sin and np.cos so that the input can be vectors
    # not just scalars.
    x1 = l1 * np.sin(theta1)
    y1 = -l1 * np.cos(theta1)
    x2 = l1 * np.sin(theta1) + l2 * np.sin(theta2)
    y2 = -l1 * np.cos(theta1) - l2 * np.cos(theta2)    
    
    return x1, y1, x2, y2



"""
Calculates the potential energy of the double pendulum.

Input:
Q: The state vector of the system. Q[0] = theta1, Q[1] = theta2, Q[2] = z1 and Q[3] = z2
M: A vector containint the mass of the bobs, M = [m1, m2]
L: A vector containing the lengths of the pedulums; L = [l1, l2]

Output: Potential energy
"""
def calculate_potential_energy(Q, M, L):
    [m1, m2] = M
    # Transform into euclidean coordinates
    [X1, Y1, X2, Y2] = polar_to_euclidean(Q, L)
    # Calculate the gravitational potential
    return m1*g*Y1 + m2*g*Y2
    
    
    
"""
Calculates the kinetic energy of the double pendulum.

Input:
Q: The state vector of the system. Q[0] = theta1, Q[1] = theta2, Q[2] = z1 and Q[3] = z2
M: A vector containint the mass of the bobs, M = [m1, m2]
L: A vector containing the lengths of the pedulums; L = [l1, l2]

Output: Kinetic energy
"""    
def calculate_kinetic_energy(Q, M, L):
    [m1, m2] = M
    [l1, l2] = L
    
    [theta1, theta2, z1, z2] = Q
     
    return 0.5*m1*l1**2*z1**2 + 0.5*m2*(l1**2*z1**2 + l2**2*z2**2 + 2*l1*l2*z1*z2*np.cos(theta1-theta2))    



"""
Calculates the total energy of the double pendulum.

Input:
Q: The state vector of the system. Q[0] = theta1, Q[1] = theta2, Q[2] = z1 and Q[3] = z2
M: A vector containint the mass of the bobs, M = [m1, m2]
L: A vector containing the lengths of the pedulums; L = [l1, l2]

Output: The total energy of the system
"""
def calculate_energy(Q, M, L):
    return calculate_kinetic_energy(Q,M,L) + calculate_potential_energy(Q,M,L)



"""
Calculates the time derivative of Q. Theis is calculated by using an analytic formula derived from
Lagranges equations.

Input:
Q: The state vector of the system. Q[0] = theta1, Q[1] = theta2, Q[2] = z1 and Q[3] = z2
M: A vector containint the mass of the bobs, M = [m1, m2]
L: A vector containing the lengths of the pedulums; L = [l1, l2]

Output: np.array containing time derivatives of theta1, theta2, z1 and z2.
"""
def derivativeVector(Q, M, L):
    [l1, l2] = L
    [m1, m2] = M
    [theta1, theta2, z1, z2] = Q
    
    sin = np.sin(theta1 - theta2)
    cos = np.cos(theta1 - theta2)
    
    
    theta1dot = z1
    theta2dot = z2
    
    z1dot = (m2*g*np.sin(theta2)*cos - m2*sin*( l1*z1**2*cos + l2*z2**2 ) - (m1 + m2)*g*np.sin(theta1) ) / (l1*( m1 + m2*sin**2))
    
    z2dot = ( (m1 + m2)*(l1*z1**2*sin - g*np.sin(theta2) + g*np.sin(theta1)*cos ) + m2*l2*z2**2*sin*cos  ) / ( l2*(m1 + m2*sin**2) )
    
    return np.array([theta1dot, theta2dot, z1dot, z2dot])
    
    
"""
Implementation of 4th order Runge-Kutta time integration. If dy/dt = f(y), then y can be 
approximated at time intervals n*dt as follows:

k1 = f(y_i)
k2 = f(y_i + dt*k1/2)
k3 = f(y_i + dt*k3/2)
k4 = f(y_i + dt*k3)

y_{i+1} = y_i + 1/6*(k1 + 2*k2 + 2*k3 + k4)

Input:
Q0: The initial configuration of the system, Q0 = [theta1_0, theta2_0, z1_0, z2_0]
M: A vector containint the mass of the bobs, M = [m1, m2]
L: A vector containing the lengths of the pedulums; L = [l1, l2]
t: An array containing all the time steps
dt: The time difference between time steps.

Output: np.array trajectory containing the configurations Q_i at each time step.

"""    
def runge_kutta(Q0,M,L,t,dt):

    start = time.time()
    
    # New values of Q are stored in trajectory.
    trajectory = np.zeros([4,len(t)])
    trajectory[:,0] = Q0
    for i in range(1,len(t)):
        # Store the configuration at the previous time step
        Q_previous = trajectory[:,i-1]
        # Calculate the intermediate steps
        k1 = derivativeVector(Q_previous, M, L)
        k2 = derivativeVector(Q_previous + dt*k1/2.0, M, L)
        k3 = derivativeVector(Q_previous + dt*k2/2.0, M, L)
        k4 = derivativeVector(Q_previous + dt*k3, M, L)
        
        # Calculate the new configuration and update trajectory
        Q_new = Q_previous + 1.0/6.0*dt*(k1 + 2*k2 + 2*k3 + k4)
        trajectory[:,i] = Q_new
    
    
    end = time.time()
    print("Runge-Kutta method running time: " + str(end-start))
    
    return trajectory
    


"""
Implementation of 4th order Adams-Bashforth method. This method is less time-consuming, since it 
doesn't calculate intermediate steps at each time step, instead it uses the derivatives at 4 earlier
steps. If dy/dt = f(y), the state can be updated as follows:

y_{i+1} = y_i + dt*(55/24*f(y_i) - 59/24*f(y_{i-1}) + 37/24*f(y_{i-2}) - 9/24*f(y_{i-3}))

Four values f(y_i) are stored in vector adams_bashforth_memory. The initial 4 entries are calculated
using Runge-Kutta.

Input:
Q0: The initial configuration of the system, Q0 = [theta1_0, theta2_0, z1_0, z2_0]
M: A vector containint the mass of the bobs, M = [m1, m2]
L: A vector containing the lengths of the pedulums; L = [l1, l2]
t: An array containing all the time steps
dt: The time difference between time steps.

Output: np.array trajectory containing the configurations Q_i at each time step.

"""
def adams_bashforth(Q0,M,L,t,dt):    

    start = time.time()

    # New values of Q are stored in trajectory.
    trajectory = np.zeros([4,len(t)])
    trajectory[:,0] = Q0
    
    # Array containing derivative vectors of 4 last configurations
    adams_bashforth_memory = np.zeros([4,4])
    # Store the derivative of the step 0
    adams_bashforth_memory[:,0] = derivativeVector(Q0, M, L)
    
    for i in range(1,4):
        # Use previous value to determine the next using Runge-Kutta 4th order method
        Q_previous = trajectory[:,i-1]
        k1 = derivativeVector(Q_previous, M, L)
        k2 = derivativeVector(Q_previous + dt*k1/2.0, M, L)
        k3 = derivativeVector(Q_previous + dt*k2/2.0, M, L)
        k4 = derivativeVector(Q_previous + dt*k3, M, L)
 
        Q_new = Q_previous + 1.0/6.0*dt*(k1 + 2*k2 + 2*k3 + k4)
        # Store the new point in trajectory
        trajectory[:,i] = Q_new    

        """
        Store the derivative at the new point in memory. In this loop the value at index i of
        adams_bashforth_memory is the derivative at the point of step i.
        """
        adams_bashforth_memory[:,i] = derivativeVector(Q_new, M, L)
    
    # Now adams_bashforth_memory = [ f(Q0), f(Q1), f(Q2), f(Q3)]
    
    """
    Now the proper Adams-Bashforth method can begin. At time step i the vector 
    adams_bashforth_memory contains values f(Q_{i-4}), f(Q_{i-3}), f(Q_{i-2}) and f(Q_1).
    The state is updated using 4th order Adams-Bashforth method:
    Q_i = Q_{i-1} + dt*( 55/24*f(Q_{i-1}) - 59/24*f(Q_{i-2}) + 37/24*f(Q_{i-3}) - 9/24*f(Q_{i-4}) )
    """    
    for i in range(4,len(t)):
        Q_previous = trajectory[:,i-1]
        """
        Here % operator is used so that the entries in adams_bashforth_memory are used cyclically. 
        This removes the need to move the entries around, thus saving processing time.
        """
        Q_new = Q_previous + dt * ( 55.0/24*adams_bashforth_memory[:,(i-1)%4] 
            - 59.0/24*adams_bashforth_memory[:,(i-2)%4] + 37.0/24*adams_bashforth_memory[:,(i-3)%4] 
            - 9.0/24*adams_bashforth_memory[:,(i-4)%4])  
        
        # Store the new state in trajectory
        trajectory[:,i] = Q_new
    
        """
        Store the derivative at the new point in memory at the correct index i modulo 4.
        """      
        adams_bashforth_memory[:,i%4] = derivativeVector(Q_new, M, L)
        
    end = time.time()
    print("Adams-Bashforth method running time: " + str(end-start))  
              
    return trajectory    
        
        
"""
Method for simulating double pendulum. 

Input:
Q0: The initial configuration of the system, Q0 = [theta1_0, theta2_0, z1_0, z2_0]
M: A vector containint the mass of the bobs, M = [m1, m2]
L: A vector containing the lengths of the pedulums; L = [l1, l2]
t: An array containing all the time steps
dt: The time difference between time steps.
method: A string specifying which method is used for the time integration. method == "RK4" or "AB4".

Output:
The output is either runge_kutta(Q0,M,L,t,dt) or adams_bashforth(Q0,M,L,t,dt) depending on the 
input parameter 'method'.

"""        
def simulate_double_pendulum(Q0,M,L,t,dt,method):
    if method == "RK4":
        return runge_kutta(Q0,M,L,t,dt)
    elif method == "AB4":
        return adams_bashforth(Q0,M,L,t,dt)
    else:
        print("Not a valid method. Valid methods are: RK4 and AB4.")
                        
        
        


"""
Method for animating the simulation of the double pendulum.

Input:
list_of_trajectories: A list containing all the trajectories to be drawn
L: Lengths of the pendulums
dt: The time difference between time steps
fps: The frames per second of the animation
"""        
def animate_double_pendulum(list_of_trajectories, L, dt, fps=30.0):
    frame_time = 1.0/fps     
    indices_in_a_frame = int(frame_time / dt)
    # The plotting area is a square with diameter twice the length of the pendulums combined.
    draw_radius = L[0] + L[1]
    
    #X1, Y1, X2, Y2 = polar_to_euclidean(trajectory, L)
    
    # Create a new figure
    fig = plt.figure()
    # Set the axes
    plt.axis([-draw_radius, draw_radius, -draw_radius, draw_radius])
    plt.gca().set_aspect('equal', adjustable='box')
    plt.axis('off')
    fig.show()
    
    # For every indices_in_a_frame:th index draw the configuration of the pendulum
    for i in range(0, len(list_of_trajectories[0][0]), indices_in_a_frame):
        fig.clear()
        
        # Loop over the list of trajectories and draw them into the same figure
        for trajectory in list_of_trajectories:
            # Transform the current trajectory into euclidean coordinates
            X1, Y1, X2, Y2 = polar_to_euclidean(trajectory, L)
            # Plot the configuration of this trajectory at time step i 
            plt.plot([0.0, X1[i], X2[i]], [ 0.0, Y1[i], Y2[i] ], '-o', ms=10)
            
        plt.axis([-draw_radius, draw_radius, -draw_radius, draw_radius])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.axis('off')
        
        fig.canvas.draw()
        #plt.savefig('frames/_img{:04d}.png'.format(i//frame_index), dpi=72)
        #plt.savefig('frames/image{:d}.png'.format(i//indices_in_a_frame), dpi=72, bbox_inches='tight')
        plt.pause(0.00001)
 
















