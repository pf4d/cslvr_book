from fenics import *

# Create mesh and define function space
mesh = UnitSquareMesh(32, 32)
V    = FunctionSpace(mesh, "Lagrange", 1)
ff   = FacetFunction('size_t', mesh, 0)

# iterate through the facets and mark each if on a boundary :
#
# 1 - West
# 2 - East
# 3 - North
# 4 - South
for f in facets(mesh):
  n       = f.normal()    # unit normal vector to facet f
  tol     = DOLFIN_EPS
  if   n.x() <= -tol and n.y() <   tol and f.exterior():
    ff[f] = 1
  elif n.x() >=  tol and n.y() <   tol and f.exterior():
    ff[f] = 2
  elif n.x() <   tol and n.y() >=  tol and f.exterior():
    ff[f] = 3
  elif n.x() <   tol and n.y() <= -tol and f.exterior():
    ff[f] = 4

ds = Measure('ds')[ff]
dN = ds(3)
dS = ds(4)
dE = ds(2)
dW = ds(1)

# Define boundary condition
u0  = Constant(1.0)
bcE = DirichletBC(V, u0, ff, 2)
bcW = DirichletBC(V, u0, ff, 1)
bc  = [bcE, bcW]

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("10 * sin(2*pi*x[0]) * sin(2*pi*x[1])")
g = Expression("sin(x[0])")
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*(dN + dS)

# Compute solution
u = Function(V)
solve(a == L, u, bc)

# get array componets and triangulation :
v = u.compute_vertex_values(mesh)
x = mesh.coordinates()[:,0]
y = mesh.coordinates()[:,1]
t = mesh.cells()

from pylab                   import *
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = figure(figsize=(8,7))
ax  = fig.add_subplot(111)

cm = get_cmap('viridis')
c  = ax.tricontourf(x, y, t, v, 10, cmap=cm)
p  = ax.triplot(x, y, t, '-', color='k', lw=0.2, alpha=0.4)
ax.axis('equal')
ax.set_xlim([x.min(), x.max()])
ax.set_ylim([y.min(), y.max()])
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_xticklabels([])
ax.set_yticklabels([])
  
divider = make_axes_locatable(gca())
cax  = divider.append_axes('right', "5%", pad="3%")
cbar = colorbar(c, cax=cax) 
tight_layout()
savefig('../../images/fenics_intro/2Dpoisson.pdf')



