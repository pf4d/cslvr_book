#===============================================================================
# finite-element solution :

from fenics import *

mesh = IntervalMesh(1000,0,1)
Q    = FunctionSpace(mesh, 'CG', 1)
x    = SpatialCoordinate(mesh)[0]
t    = mesh.coordinates()[:,0]
n    = FacetNormal(mesh)

v    = TestFunction(Q)
u    = Function(Q)
du   = TrialFunction(Q)

eps  = 0.05

dur  = -10.0

r    = - eps * u.dx(0) * v.dx(0) * dx \
       + eps * dur * v * ds \
       - (2*x + 1) * u.dx(0) * v * dx  \
       + 2 * u * v * dx \

def left(x, on_boundary):
  return on_boundary and x[0] == 0.0

def right(x, on_boundary):
  return on_boundary and x[0] == 1.0

leftBC  = DirichletBC(Q, 1.0, left)
rightBC = DirichletBC(Q, 0.0, right)

bcs = [leftBC]
J   = derivative(r, u, du)

solve(r == 0, u, bcs=bcs, J=J)

uf = u.vector().array()[::-1]

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

sf  = dur*(t - t[-1]) + uf[-1]

ax.set_xlim(0.0, 1.00)
ax.set_ylim(0.0, 3.05)

ax.plot(t, uf, 'k',   lw=2.0, label=r"FEniCS")
ax.plot(t, sf, 'r--', lw=2.0, label=r"slope")

axins = zoomed_inset_axes(ax, 8, loc=4)
axins.set_xlim(0.95, 1.00)
axins.set_ylim(2.75, 3.00)
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
xticks(visible=False)
yticks(visible=False)

axins.plot(t, uf, 'k',   lw=2.0, label=r"FEniCS")
axins.plot(t, sf, 'r--', lw=2.0, label=r"slope")
axins.grid()

leg = ax.legend(loc='upper left', fontsize='medium')
leg.get_frame().set_alpha(0.0)
ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$u$')
ax.grid()

tight_layout()
savefig("../../images/fenics_intro/1D_BVP_2.pdf")



