from fenics import *

mesh = IntervalMesh(1000,0,2*pi)
Q    = FunctionSpace(mesh, 'CG', 1)
x    = SpatialCoordinate(mesh)[0]
t    = mesh.coordinates()[:,0]

u    = interpolate(Expression('cos(x[0])'), Q)
s    = interpolate(Expression('sin(x[0])'), Q)
v    = TrialFunction(Q)
phi  = TestFunction(Q)

def left(x, on_boundary):
  return on_boundary and abs(x[0]) < 1e-14

# integral is zero on the left
bcs = DirichletBC(Q, 0.0, left)

a      = v.dx(0) * phi * dx
L      = u * phi * dx
v      = Function(Q)
solve(a == L, v, bcs)

uf = u.vector().array()[::-1]
vf = v.vector().array()[::-1]
sf = s.vector().array()[::-1]

r  = vf - sf

#===============================================================================
# plot :

from pylab import *
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

mpl.rcParams['font.family']          = 'serif'
mpl.rcParams['text.usetex']          = True
mpl.rcParams['text.latex.preamble']  = ['\usepackage{fouriernc}']

fig = figure(figsize=(5,3.5))
ax  = fig.add_subplot(111)

ax.plot(t, uf, 'k',       lw=2.0, label=r"$u$")
ax.plot(t, vf, 'r--',     lw=2.0, label=r"$v$")
ax.plot(t, r,  '#880cbc', lw=2.0, label=r"$\epsilon$")

ax.set_xlim(0,2*pi)
leg = ax.legend(loc='upper center', fontsize='medium')
leg.get_frame().set_alpha(0.0)
ax.set_xlabel(r'$x$')
ax.grid()

tight_layout()
savefig("../../images/fenics_intro/1D_BVP_4.pdf")


