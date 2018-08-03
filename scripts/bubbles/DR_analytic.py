from fenics import *

mesh = IntervalMesh(10,0,1)

Q    = FunctionSpace(mesh, 'CG', 1)
B    = FunctionSpace(mesh, 'B',  2)
M    = Q + B

def left(x, on_boundary):
  return on_boundary and x[0] == 0

kappa = Constant(1.0/500.0)
s     = Constant(1.0)
f     = Constant(0.0)

#==============================================================================
# standard Galerkin solution :

leftBC  = DirichletBC(Q, 1.0, left)

u   = TrialFunction(Q)
v   = TestFunction(Q)
us  = Function(Q)
uf1 = Function(Q)

a  = + kappa * u.dx(0) * v.dx(0) * dx \
     + s * u * v * dx
L  = f * v * dx

solve(a == L, us, leftBC)

uf1.interpolate(us)

#==============================================================================
# enriched space solution :

uD      = project(Constant(1.0), M)
leftBC  = DirichletBC(M, uD, left)

u   = TrialFunction(M)
v   = TestFunction(M)
us  = Function(M)
uf2 = Function(Q)

a  = + kappa * u.dx(0) * v.dx(0) * dx \
     + s * u * v * dx
L  = f * v * dx

solve(a == L, us, leftBC)

uf2.interpolate(us)

#==============================================================================
# SSM stabilized :

leftBC  = DirichletBC(Q, 1.0, left)

u   = TrialFunction(Q)
v   = TestFunction(Q)
us  = Function(Q)
uf3 = Function(Q)
h   = CellSize(mesh)
C   = Constant(1/15.0)
tau = C * h**2 / kappa

def L(u):  return -(kappa * u.dx(0)).dx(0) + s*u 

a  = + kappa * u.dx(0) * v.dx(0) * dx \
     + s * u * v * dx \
     - inner(L(v), tau*L(u)) * dx
L  = + f * v * dx \
     - inner(L(v), tau*f) * dx

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

x    = linspace(0, 1, 1000)
ue   = (exp(-10*sqrt(5)*(x-2)) + exp(10*sqrt(5)*x))/(1 + exp(20*sqrt(5)))

mpl.rcParams['font.family']          = 'serif'
mpl.rcParams['text.usetex']          = True
mpl.rcParams['text.latex.preamble']  = ['\usepackage{fouriernc}']

fig = figure(figsize=(5,3.5))
ax   = fig.add_subplot(111)

ax.plot(t, uf1, 'r',   ls='-',  lw=2.0, label=r"$u$")
ax.plot(t, uf2, 'k',   ls='-',  lw=2.0, label=r"$\hat{u}$")
ax.plot(t, uf3, purp,  ls='--', lw=2.0, label=r"$\tilde{u}$")
ax.plot(x, ue,  'k--', ls='--', lw=2.0, label=r"$u_{\mathrm{a}}$")

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$u$')
ax.grid()
leg = ax.legend(loc='upper right')
leg.get_frame().set_alpha(0.0)
ax.set_xlim([0,0.3])
tight_layout()
savefig('../../images/bubbles/DR_analytic_new.pdf')



