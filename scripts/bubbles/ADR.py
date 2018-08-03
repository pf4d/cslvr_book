from fenics import *

mesh = IntervalMesh(50,0,1)

Q    = FunctionSpace(mesh, 'CG', 1)
B    = FunctionSpace(mesh, 'B',  2)
M    = Q + B

def left(x, on_boundary):
  return on_boundary and x[0] == 0

uD        = project(Constant(0.0), M)
leftBC    = DirichletBC(Q, 0.0, left)
leftBC_b  = DirichletBC(Q, uD,  left)

kappa   = Constant(1.0/100.0)
s       = Constant(5.0)
d       = Constant(10.0)
f       = Function(Q)
f.vector()[25] = 1000  # this is about the middle for a 50 element mesh

#==============================================================================
# standard Galerkin solution :

u   = TrialFunction(Q)
v   = TestFunction(Q)
us  = Function(Q)
uf1 = Function(Q)

a  = + kappa * u.dx(0) * v.dx(0) * dx \
     + d * u.dx(0) * v * dx \
     + s * u * v * dx
L  = f * v * dx

solve(a == L, us, leftBC)

uf1.interpolate(us)

#==============================================================================
# bubble-enriched solution :
u   = TrialFunction(M)
v   = TestFunction(M)
us  = Function(M)
uf2 = Function(Q)

a  = + kappa * u.dx(0) * v.dx(0) * dx \
     + d * u.dx(0) * v * dx \
     + s * u * v * dx
L  = f * v * dx

solve(a == L, us, leftBC_b)

uf2.interpolate(us)

#==============================================================================
# GLS stabilized solution :
u   = TrialFunction(Q)
v   = TestFunction(Q)
us  = Function(Q)
uf3 = Function(Q)
h   = CellSize(mesh)

# for SUPG :
#Pe  = h * d / (2*kappa)
#tau = h / (2*d) * (1/tanh(Pe) - 1 / Pe)

# for GLS or SSM :
tau = 1 / (4*kappa/h**2 + 2*d/h + s)

def L(u):       return -(kappa * u.dx(0)).dx(0) + d*u.dx(0) + s*u  # GLS
def L_star(u):  return -(kappa * u.dx(0)).dx(0) - d*u.dx(0) + s*u  # SSM
def L_adv(u):   return d*u.dx(0)                                   # SUPG


a  = + kappa * u.dx(0) * v.dx(0) * dx \
     + d * u.dx(0) * v * dx \
     + s * u * v * dx \
     + inner(L(v), tau*L(u)) * dx
L  = + f * v * dx \
     + inner(L(v), tau*f) * dx

solve(a == L, us, leftBC)

uf3.interpolate(us)

#==============================================================================
# plotting :
from pylab import *

purp = '#880cbc'

t    = mesh.coordinates()[:,0][::-1]
uf1  = uf1.vector().array()
uf2  = uf2.vector().array()
uf3  = uf3.vector().array()

mpl.rcParams['font.family']          = 'serif'
mpl.rcParams['text.usetex']          = True
mpl.rcParams['text.latex.preamble']  = ['\usepackage{fouriernc}']

fig = figure(figsize=(5,3.5))
ax   = fig.add_subplot(111)

ax.plot(t, uf1, 'r',   ls='-', lw=2.0, label=r"$u$")
ax.plot(t, uf2, 'k',   ls='-', lw=2.0, label=r"$\hat{u}$")
ax.plot(t, uf3,  purp, ls='-', lw=2.0, label=r"$\tilde{u}$")

ax.set_xlabel(r'$x$')
ax.grid()
leg = ax.legend(loc='upper left')
leg.get_frame().set_alpha(0.0)
tight_layout()
savefig('../../images/bubbles/extreme_new.pdf')



