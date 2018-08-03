from fenics import *

mesh = IntervalMesh(1000,0,10)
Q    = FunctionSpace(mesh, 'CG', 1)

v    = interpolate(Expression('10*exp(-0.5*pow(x[0] - 5.0, 2)/pow(2,2))'), Q)
u    = 1.0
g    = v / u

y    = Function(Q)
dy   = TrialFunction(Q)
phi  = TestFunction(Q)

def left(x, on_boundary):
  return on_boundary and abs(x[0]) < 1e-14

# boat starts at y=0 :
bc = DirichletBC(Q, 0.0, left)

def Lp(yp):
  Lp  = yp * (1 + yp**2)**(-1/2.) / sqrt(v**2 + u**2)
  return Lp

def Lpp(yp):
  Lpp = 1/(sqrt(1 + yp**2) * sqrt(v**2 + u**2)) * (1 - yp**2 / (1 + yp**2)**3)
  return Lpp

def M(yp):
  q = (2*yp**2 + 1) / (sqrt(u**2 + v**2) * sqrt(yp**2 + 1))
  return q

F  = + (Lpp(y.dx(0)) + M(y.dx(0))) * y.dx(0) * phi.dx(0) * dx \
     - (M(g) + Lpp(g))* g * phi * ds \
     + Lp(g) * phi * ds

J  = derivative(F, y, dy)

solve(F == 0, y, bc, J=J)

#===============================================================================
# plot :

from pylab import *

mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['legend.fontsize'] = 'medium'

x  = mesh.coordinates()[:,0][::-1]
vf = v.vector().array()
yf = y.vector().array()

fig = figure(figsize=(6,4))
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()

ax1.plot(x, yf, 'k',   lw=2.0, label=r"$y$")
ax2.plot(x, vf, 'r--', lw=2.0, label=r"$v$")

ax2.tick_params(axis='x', colors='r')
ax2.yaxis.label.set_color('r')
ax2.set_ylabel(r'$v$')
leg2 = ax2.legend(loc='upper right')

leg2.get_frame().set_alpha(0.0)

for tl in ax2.get_yticklabels():
  tl.set_color('r')

#ax1.set_ylim(0,2*pi)
ax1.set_ylim(0,10)
ax1.set_xlabel(r'$x$')
ax1.set_ylabel(r'$y$')
ax1.grid()
leg1 = ax1.legend(loc='upper left')
leg1.get_frame().set_alpha(0.0)

tight_layout()
savefig("../../images/fenics_intro/river_cross.pdf")
show()
