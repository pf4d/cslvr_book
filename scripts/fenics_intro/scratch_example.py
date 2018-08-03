from pylab        import *
from fenics       import *
from scipy.linalg import lu
from scipy.linalg import solve as sl

#===============================================================================
# manual solution :

K  = array([[6,-3,0],[-3,6,-3],[0,-3,3]])
F  = array([1/3., 1/3., 1/6.])

P,L,U = lu(K)
y     = sl(L,F)
uf    = sl(U,y)
uf    = append(0.0, uf)

#===============================================================================
# FEniCS solution :

mesh = IntervalMesh(3,0,1)
Q    = FunctionSpace(mesh, "CG", 1)

f    = Constant(1.0)
u    = TrialFunction(Q)
v    = TestFunction(Q)

a    = u.dx(0) * v.dx(0) * dx
l    = f * v * dx

def left(x, on_boundary):
  return x[0] == 0 and on_boundary
bc   = DirichletBC(Q, 0.0, left)

u    = Function(Q)
solve(a == l, u, bc)

#===============================================================================
# plotting :

xe = linspace(0,1,1000)
xf = linspace(0,1,4)

uv = u.vector().array()[::-1]
ue = -0.5*(xe - 2)*xe
us = xe

mpl.rcParams['font.family']          = 'serif'
mpl.rcParams['text.usetex']          = True
mpl.rcParams['text.latex.preamble']  = ['\usepackage{fouriernc}']

fig = figure(figsize=(5,3))
ax  = fig.add_subplot(111)

ax.plot(xe, ue, 'k--', lw=2.0, label='exact solution')
ax.plot(xf, uv, 'k-',  lw=2.0, label='FEniCS solution')
ax.plot(xf, uf, 'ro',  lw=2.0, label='solution by hand')
ax.plot(xe, us, 'r--', lw=2.0, label='flux of $u$ on left')

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$u$')
ax.set_ylim([0,0.6])
ax.grid()
leg = ax.legend(loc='lower right', fontsize='medium')
leg.get_frame().set_alpha(0.0)
tight_layout()
savefig("../../images/fenics_intro/scratch_example.pdf")



