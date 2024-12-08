# This is a copy implementation from the third video from Paul's physics 
# lectures on YouTube. This needs to be modified for the 2D case.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# Parameters
L       = 10.0      # length of "string"
x0      = 5.0       # Starting point of wave
y0      = 5.0       # Starting point of wave
k       = 10        # parameter for I(x) eq
Tf      = 25.0      # Time duration
v       = 1.0       # Wave velocity
n       = 1000      # Number of grid points
t_steps = 400  # Number of time steps

# Grid for graph of displacement versus time
x, dx   = np.linspace(0, L, n, retstep = True)
y, dy   = np.linspace(0, L, n, retstep = True)
x,y = np.meshgrid(x,y)
t, dt   = np.linspace(0, Tf, t_steps, retstep = True)
print(dx, dy, dt)

# Assume step sizes for x and y are the same
# i.e. dx = dy
r   = (v*dt)/dx     # Courant number; used in finite diff method
c   = r**2          # usefule quantity

print(r)

# Gaussian pluck; 
# this equation defines the boundary conditions for all x at t=0
def Io(x,y):
    #u = np.zeros((len(x), len(y)))
    #for i in range(len(x)):
    #    for j in range(len(y)):
    #        u[i,j] = np.exp(-k*((x[i]-x0)**2+(y[j]-y0)**2))

    u = np.exp(-k*((x-x0)**2+(y-y0)**2))
    return u

# du/dx boundary condition
# at t=0, du/dx=0 and du/dy=0 for all (x,y)
def g(x,y):
    G = 0.0
    return G

# array for storing U(x,y;t) 
U = np.zeros((len(x), len(y), len(t)))

# Assign boundary conditions in U(x,t)
U[0,:, :]     = 0.0   # U(0,y;t) = 0.0
U[-1,:, :]    = 0.0   # U(L,y;t) = 0.0
U[:,0, :]     = 0.0   # U(x,0;t) = 0.0
U[:,-1, :]    = 0.0   # U(x,L;t) = 0.0
U[:,:, 0]     = Io(x,y) # at t=0, U(x,y;0)=Io(x,y)

# Using derived finite difference method formulas

print("Calculate initial values ...")
# define values for U(x_i, 1) i.e. values just after boundaries
for XX in range(1, len(x)-1):
    for YY in range(1, len(x)-1):
        U[XX, YY, 1] = 0.5*(c*U[XX+1,YY,0] + c*U[XX-1,YY,0] + c*U[XX,YY+1,0] + c*U[XX,YY-1,0] + 2*(1-2*c)*U[XX,YY,0]) + dt*g(x,y)
        #print(XX, YY, U[XX, YY, 1])

print("Calculate successive values ...")
# define values for U(x_i, t) for all values t>1
for tt in range(1, len(t)-1):
    for XX in range(1, len(x)-1):
        for YY in range(1, len(x)-1):
            U[XX, YY, tt+1] = (c*U[XX+1,YY,tt] + c*U[XX-1,YY,tt] + c*U[XX,YY+1,tt] + c*U[XX,YY-1,tt] + 2*(1-2*c)*U[XX,YY,tt]) -U[XX, YY, tt-1]
        #print(XX, YY, U[XX, YY, tt+1])



print("Time to plot!!")

#x,y = np.meshgrid(x,y)

# plot results
fig = plt.figure()
ax1 = fig.add_subplot(2,1,1, projection='3d')

show_time = 50 

ax1.plot_surface(x,y,U[:,:,show_time], rstride=3, cstride=3, linewidth=1, antialiased=True, cmap=cm.viridis)
ax1.view_init(55,-70)
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_zticks([])
ax1.set_xlabel("x axis")
ax1.set_ylabel("y axis")

ax2 = fig.add_subplot(2,1,2, projection='3d')

ax2.contourf(x,y,U[:,:,show_time], zdir='z', offset=0, cmap=cm.viridis)
ax2.grid(False)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_zticks([])
ax2.set_xlabel("x axis")
ax2.set_ylabel("y axis")


plt.show()
