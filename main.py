'''
rho = [0,0,0]
rho_s = [0,0,0]

sp_dot = 0

Only translational dynamics are considered
'''
import numpy as np
from RK4_project import RK4 as RK4

# define simulation time
start = 0
end = 10
N_time = 100*(end-start)+1
time = np.linspace(start, end, N_time)
h = time[1]-time[0]

g = 0.0 # gravity
L = 1.5 # length of the cable
EA = 659700.0 # value of EA, copied from the reference
mb = 0.7 # total mass on the quadrotor in kg
ms = 0.5 # total mass of the payload

# define initial states
x0 = np.zeros(3) 
dxdt0 = np.array([0.0,0.0,1.0]) # dxdt0

rL0 = np.array([0.0,0.0,L]) # rL0
drLdt0 = np.array([0.0,0.0,1.0]) # drLdt0

# define a clubbed state
X_array = np.zeros((12, N_time))

# initialization
X_array[0:3, 0] = x0
X_array[3:6, 0] = dxdt0
X_array[6:9, 0] = rL0
X_array[9:12, 0] = drLdt0

def getStates(X_array):
    x       = X_array[0:3]
    dxdt    = X_array[3:6]
    rL      = X_array[6:9]
    drLdt   = X_array[9:12]
    return x, dxdt, rL, drLdt

for i in range(N_time-1):
    print(X_array[:,i])

    # define the state space equation X_dot = f(t, X) 
    def dynamics(t, X_array):
        x       = getStates(X_array)[0]
        dxdt    = getStates(X_array)[1]
        rL      = getStates(X_array)[2]
        drLdt   = getStates(X_array)[3]
        r_prime = (rL - x)/L
        norm    = np.linalg.norm(r_prime)

        F = EA*(norm-1.0)*r_prime/norm # tension force

        X_dot = np.zeros(12)
        X_dot[0:3] = dxdt
        X_dot[3:6] = g*np.array([0.0,0.0,1.0]) + F/mb
        X_dot[6:9] = drLdt
        X_dot[9:12] = g*np.array([0.0,0.0,1.0]) - F/ms
        return X_dot
    
    X_dot = dynamics(0, X_array[:,i])
    
    # integration using RK4
    X_array[:,i+1] = RK4(dynamics, X_array[:,i], time[i], time[i+1], h/10.0)
    
    #integration using Euler's method
    # X_array[:,i+1] = X_array[:,i] + h*X_dot