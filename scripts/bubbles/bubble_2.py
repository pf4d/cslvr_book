from pylab                 import *
from mpl_toolkits.mplot3d  import Axes3D

def B(x,y):
  return 27 * (1-x-y) * x * y

x   = linspace(0,1,100)
y   = linspace(0,1,100)

X,Y = meshgrid(x,y)

# plot the results :
fig = plt.figure()
ax  = fig.add_subplot(111, projection='3d')

Z   = B(X,Y)

Z[X + Y > 1] = 0  # zero the area outside the triangle

ax.plot_wireframe(X, Y, Z, color='r', lw=2.0, rstride=5, cstride=5)

ax.set_ylabel(r'$y$')
ax.set_xlabel(r'$x$')
ax.set_zlabel(r'$B$')
tight_layout()
show() 
