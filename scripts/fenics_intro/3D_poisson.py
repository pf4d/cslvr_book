from fenics import *

# create mesh :
p1   = Point(-pi, -pi, -pi)
p2   = Point( pi,  pi,  pi)
mesh = BoxMesh(p1, p2, 30, 30, 30)

# define function space :
V    = FunctionSpace(mesh, "Lagrange", 1)

# create a MeshFunction for marking boundaries :
ff   = FacetFunction('size_t', mesh, 0)

# iterate through the facets and mark each if on a boundary :
#
#   1 =  ..... top           |       4 =  ..... West side
#   2 =  ..... bottom        |       5 =  ..... North side
#   3 =  ..... East side     |       6 =  ..... South side
for f in facets(mesh):
  n       = f.normal()    # unit normal vector to facet f
  tol     = 1e-10
  if   n.z() >=  tol and f.exterior():
    ff[f] = 1
  elif n.z() <= -tol and f.exterior():
    ff[f] = 2
  elif n.x() > tol   and n.y() < tol and f.exterior():
    ff[f] = 3
  elif n.x() < -tol  and n.y() < tol and f.exterior():
    ff[f] = 4
  elif n.y() > tol   and n.x() < tol and f.exterior():
    ff[f] = 5
  elif n.y() < -tol  and n.x() < tol and f.exterior():
    ff[f] = 6

# need the N and S boundary for natural conditions :
ds   = Measure('ds')[ff]
dN   = ds(5) + ds(6)

# define essential boundary conditions :
zero = Constant(0.0)
bcN  = DirichletBC(V, zero, ff, 2)
bcS  = DirichletBC(V, zero, ff, 1)

# define variational functions :
u = TrialFunction(V)
v = TestFunction(V)

# expressions for known data :
f = Expression("10 * exp(-(pow(x[0],2)/2 + pow(x[1],2)/2 + pow(x[2],2)/2))")
g = Expression("sin(x[0])")

# variational problem :
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*dN

# compute solution :
u = Function(V)
solve(a == L, u, [bcN, bcS])

File("output/u.pvd") << u


