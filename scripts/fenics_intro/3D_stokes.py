from fenics import *

mesh    = Mesh('meshes/unit_cyl_mesh.xml')

# Define function spaces
#B  = FunctionSpace(mesh, "B", 4)
#Q  = FunctionSpace(mesh, "CG", 1)
#M  = Q + B
#V  = MixedFunctionSpace([M,M,M])
#W  = MixedFunctionSpace([V,Q])
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)
W = V * Q
ff = FacetFunction('size_t', mesh, 0)

# iterate through the facets and mark each if on a boundary :
#
#   1 = high slope, upward facing ................ top
#   2 = high slope, downward facing .............. bottom
#   3 = low slope, upward or downward facing ..... side
for f in facets(mesh):
  n       = f.normal()    # unit normal vector to facet f
  tol     = 1.0
  if   n.z() >=  tol and f.exterior():
    ff[f] = 1
  elif n.z() <= -tol and f.exterior():
    ff[f] = 2
  elif abs(n.z()) < tol and f.exterior():
    ff[f] = 3

L    = 1.0/2000.0
xmin = -L
xmax = L
ymin = -L
ymax = L

# Deform the mesh to the defined geometry :
for x in mesh.coordinates():
  # transform x :
  x[0]  = x[0]  * (xmax - xmin)
  # transform y :
  x[1]  = x[1]  * (ymax - ymin)
  # transform z :
  # thickness = surface - base, z = thickness + base
  x[2]  = x[2] * 5.0/1000.0

# constants :
rho    = Constant(1420.0)
eta    = Constant(8.0)
g      = Constant(9.8)
x      = SpatialCoordinate(mesh)
n      = FacetNormal(mesh)
I      = Identity(3)

#===============================================================================
# define variational problem :
U    = TrialFunction(W)
tst  = TestFunction(W)

u, p = split(U)
v, q = split(tst)

# no-slip boundary condition for velocity :
bc   = DirichletBC(W.sub(0), Constant((0, 0, 0)), ff, 3)

# stress and strain tensors :
def sigma(u,p): return 2*eta*epsilon(u) - p*I
def epsilon(u): return 0.5*(grad(u) + grad(u).T)

# internal force vector :
f   = rho * as_vector([0, 0, -g])

# conservation of momentum :
a = + inner(sigma(u,p), grad(v)) * dx \
    + div(u) * q * dx
L = + dot(f, v) * dx

# solve the linear system :
U = Function(W)
solve(a == L, U, bc, solver_parameters = {"linear_solver"  : "minres", 
                                          "preconditioner" : "hypre_amg"} )
u,p = U.split(True)

File("output/u.pvd") << u
File("output/p.pvd") << p



