# This is a copy implementation from the third video from Paul's physics 
# lectures on YouTube. This needs to be modified for the 2D case.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


# The main parameters that contribute to execution time
n       = 100      # Number of grid points
t_steps = 200  # Number of time steps


# Parameters
L       = 10.0      # length of "string"
x0      = 5.0       # Starting point of wave
y0      = 5.0       # Starting point of wave
k       = 10        # parameter for I(x) eq
Tf      = 10.0      # Time duration
v       = 1.0       # Wave velocity

# Grid for graph of displacement versus time
x, dx   = np.linspace(0, L, n, retstep = True)
y, dy   = np.linspace(0, L, n, retstep = True)
np.savetxt("x_coords_premesh.csv", x, delimiter=',')

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


# Output data to text file
np.savetxt("u_coords.csv", U[:,:,0], delimiter=',')
np.savetxt("x_coords.csv", x, delimiter=',')
np.savetxt("y_coords.csv", y, delimiter=',')

f = open("x_coords.txt", "a")
f.write("Step size = " + str(dx))
f.close()

f = open("y_coords.txt", "a")
f.write("Step size = " + str(dy))
f.close()

print("Time to plot!!")

#x,y = np.meshgrid(x,y)

# plot results
fig = plt.figure()
ax1 = fig.add_subplot(2,1,1, projection='3d')
#ax2 = fig.add_subplot(2,2,1, projection='3d')
#ax3 = fig.add_subplot(2,1,2, projection='3d')

show_time = 20 

for T in range(0,len(t), 10):

    ax1.plot_surface(x,y,U[:,:,T], rstride=3, cstride=3, linewidth=1, antialiased=True, cmap=cm.viridis)

ax1.view_init(55,-70)
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_zticks([])
ax1.set_xlabel("x axis")
ax1.set_ylabel("y axis")

#for T in range(0,50,10):
#
#    ax2.plot_surface(x,y,U[:,:,T], rstride=3, cstride=3, linewidth=1, antialiased=True, cmap=cm.viridis)
#
#ax2.view_init(55,-70)
#ax2.set_xticks([])
#ax2.set_yticks([])
#ax2.set_zticks([])
#ax2.set_xlabel("x axis")
#ax2.set_ylabel("y axis")
#
#
#ax3.plot_surface(x,y,U[:,:,1], rstride=3, cstride=3, linewidth=1, antialiased=True, cmap=cm.viridis)
#for T in range(50,100,10):
#
#    ax3.plot_surface(x,y,U[:,:,T], rstride=3, cstride=3, linewidth=1, antialiased=True, cmap=cm.viridis)
#
#ax3.view_init(55,-70)
#ax3.set_xticks([])
#ax3.set_yticks([])
#ax3.set_zticks([])
#ax3.set_xlabel("x axis")
#ax3.set_ylabel("y axis")


#ax1_c = fig.add_subplot(2,1,1, projection='3d')
#ax2_c = fig.add_subplot(2,1,2, projection='3d')
#ax3_c = fig.add_subplot(2,2,1, projection='3d')

#ax1_c.contourf(x,y,U[:,:,0], zdir='z', offset=0, cmap=cm.viridis)
#
#ax1_c.grid(False)
#ax1_c.set_xticks([])
#ax1_c.set_yticks([])
#ax1_c.set_zticks([])
#ax1_c.set_xlabel("x axis")
#ax1_c.set_ylabel("y axis")
#
#ax2_c.contourf(x,y,U[:,:,50], zdir='z', offset=0, cmap=cm.viridis)
#
#ax2_c.grid(False)
#ax2_c.set_xticks([])
#ax2_c.set_yticks([])
#ax2_c.set_zticks([])
#ax2_c.set_xlabel("x axis")
#ax2_c.set_ylabel("y axis")
#
#
#ax3_c.contourf(x,y,U[:,:,200], zdir='z', offset=0, cmap=cm.viridis)
#
#ax3_c.grid(False)
#ax3_c.set_xticks([])
#ax3_c.set_yticks([])
#ax3_c.set_zticks([])
#ax3_c.set_xlabel("x axis")
#ax3_c.set_ylabel("y axis")

plt.show()
