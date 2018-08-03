#===============================================================================
# finite-element solution :

from fenics import *

mesh = IntervalMesh(1000,0,1)
Q    = FunctionSpace(mesh, 'CG', 1)
x    = SpatialCoordinate(mesh)[0]
t    = mesh.coordinates()[:,0]

v    = TestFunction(Q)
u    = Function(Q)
du   = TrialFunction(Q)

eps  = 0.05

r    = - eps * u.dx(0) * v.dx(0) * dx \
       - (2*x + 1) * u.dx(0) * v * dx  \
       + 2 * u * v * dx \

def left(x, on_boundary):
  return on_boundary and x[0] == 0.0

def right(x, on_boundary):
  return on_boundary and x[0] == 1.0

leftBC  = DirichletBC(Q, 1.0, left)
rightBC = DirichletBC(Q, 0.0, right)

bcs = [leftBC, rightBC]
J   = derivative(r, u, du)

solve(r == 0, u, bcs=bcs, J=J)

uf = u.vector().array()[::-1]

#===============================================================================
# singular-perturbation method solution :

from pylab  import *

uo = 2*t + 1
ui = -3*(exp((3*t - 3)/eps) - 1)
ul = 3

uu = uo + ui - ul

#===============================================================================
# plot :

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

mpl.rcParams['font.family']          = 'serif'
mpl.rcParams['text.usetex']          = True
mpl.rcParams['text.latex.preamble']  = ['\usepackage{fouriernc}']

fig = figure(figsize=(5,3.5))
ax  = fig.add_subplot(111)

ax.set_xlim(0.0, 1.05)
ax.set_ylim(0.0, 3.05)

ax.plot(t, uf, 'k',   lw=2.0, label=r"FEniCS")
ax.plot(t, uo, 'r--', lw=2.0, label=r"$u_o(t)$")
ax.plot(t, ui, 'r:',  lw=2.0, label=r"$u_o(t)$")
ax.plot(t, uu, 'k--', lw=2.0, label=r"$u_u(t)$")

axins = zoomed_inset_axes(ax, 4, loc=8)
axins.set_xlim(0.88, 0.98)
axins.set_ylim(2.60, 2.95)
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
xticks(visible=False)
yticks(visible=False)

axins.plot(t, uf, 'k',   lw=2.0, label=r"FEniCS")
axins.plot(t, uo, 'r--', lw=2.0, label=r"$u_o(t)$")
axins.plot(t, ui, 'r:',  lw=2.0, label=r"$u_o(t)$")
axins.plot(t, uu, 'k--', lw=2.0, label=r"$u_u(t)$")
axins.grid()

leg = ax.legend(loc='upper left', ncol=2, fontsize='medium')
leg.get_frame().set_alpha(0.0)
ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$u$')
ax.grid()

tight_layout()
savefig("../../images/fenics_intro/1D_BVP_1.pdf")




