# This is a copy implementation from the third video from Paul's physics 
# lectures on YouTube. This needs to be modified for the 2D case.

import numpy as np
import matplotlib.pyplot as plt

# Parameters
L   = 10.0      # length of "string"
x0  = 5.0       # Starting point of wave
k   = 10        # parameter for I(x) eq
Tf  = 25.0      # Time duration
v   = 1.0       # Wave velocity

# Grid for graph of displacement versus time
x, dx   = np.linspace(0, L, 1001, retstep = True)
t, dt   = np.linspace(0, Tf, 4001, retstep = True)
print(dx, dt)

r   = (v*dt)/dx     # Courant number; used in finite diff method
c   = r**2          # usefule quantity

print(r)

# Gaussian pluck; 
# this equation defines the boundary conditions for all x at t=0
def Io(x):
    u = np.exp(-k*(x-x0)**2)
    return u

# du/dx boundary condition
# at t=0, du/dx=0 for all x
def g(x):
    G = 0.0
    return G

# array for storing U(x,t) 
U = np.zeros((len(x), len(t)))

# Assign boundary conditions in U(x,t)
U[0, :]     = 0.0   # U(0,t) = 0.0
U[-1, :]    = 0.0   # U(L,t) = 0.0
U[:, 0]     = Io(x) # at t=0, U(x)=Io(x)

# Using derived finite difference method formulas

print("Calculate initial values ...")
# define values for U(x_i, 1) i.e. values just after boundaries
for XX in range(1, len(x)-1):
    U[XX, 1] = 0.5*(c*U[XX+1,0] + 2*(1-c)*U[XX,0] + c*U[XX-1, 0]) + dt*g(x)

print("Calculate successive values ...")
# define values for U(x_i, t) for all values t>1
for tt in range(1, len(t)-1):
    for XX in range(1, len(x)-1):
        U[XX, tt+1] = c*U[XX+1,tt] + 2*(1-c)*U[XX,tt] + c*U[XX-1, tt] - U[XX, tt-1] 


print("Time to plot!!")

# plot results
for T in range(0, len(t), 300):
    plt.plot(x, U[:,T], ls = '-.', lw = 1, label = 't = ' +str(T))

plt.legend()
plt.tight_layout()
plt.grid()

plt.show()
