from pylab  import *
from fenics import *
import sympy as sp

n    = 8                                 # number of vertices in mesh
l    = 1                                 # length of the domain
mesh = IntervalMesh(n, 0, l)             # 1D mesh
Q    = FunctionSpace(mesh, 'CG', 1)      # linear-Lagrange F.S.
phi  = TestFunction(Q)
u    = TrialFunction(Q)
tol  = 1e-14

def left(x, on_boundary):  return on_boundary and abs(x[0]) < tol

def right(x, on_boundary): return on_boundary and abs(x[0] - l) < tol

gamma_l = DirichletBC(Q, 0.0, left)      # left boundary condition
gamma_r = DirichletBC(Q, 0.0, right)     # right boundary condition

k = inner(grad(phi), grad(u)) * dx       # weak 2nd derivative matrix
m = phi * u * dx                         # mass matrix

K = PETScMatrix()                        # container for stiffness matrix
M = PETScMatrix()                        # container for mass matrix

K = assemble(k, tensor=K)                # global assembly of stiffness mat.
M = assemble(m, tensor=M)                # global assembly of mass mat.

gamma_r.apply(K)                         # apply right b.c. to K
gamma_l.apply(K)                         # apply left b.c. to K
gamma_r.apply(M)                         # apply right b.c. to M
gamma_l.apply(M)                         # apply left b.c. to M

eigensolver = SLEPcEigenSolver(K,M)      # create solver for Kx = lambda Mx
eigensolver.parameters['solver'] = 'lapack'
eigensolver.solve()                      # solve with LAPACK

# generate the results and calculate the Fourier coef's with SymPy
x     = sp.symbols('x')
xm    = mesh.coordinates()[:,0][::-1]
xf    = linspace(0,l,1000)
col   = ['k', 'r', '#880cbc']
lam_m = []
lam_f = []
Am    = []
Af    = []
Av    = []
cn    = []

for i in range(n-1):
  ep     = eigensolver.get_eigenpair(n-2-i)
  Afi    = zeros(n)
  lam_fi = ((i+1) * pi)**2
  Afi    = np.sin(np.sqrt(lam_fi) * xf)
  cni    = 2/l * sp.integrate(sp.sin((i+1)*pi*x/l), (x,0,l))
  lam_m.append(ep[0])
  lam_f.append(lam_fi)
  e_v = ep[2].array()
  if i == 1:
    e_v = -e_v
  Am.append(e_v)
  Af.append(Afi)
  cn.append(cni)

# plotting :
#===============================================================================

mpl.rcParams['font.family']          = 'serif'
mpl.rcParams['text.usetex']          = True
mpl.rcParams['text.latex.preamble']  = ['\usepackage{fouriernc}']

fig = figure(figsize=(5,3.5))
ax  = fig.add_subplot(111)

for i, (am,af,lf,lm,c) in enumerate(zip(Am, Af, lam_f, lam_m, col)):
  lbl = r'$n=%i,\ \lambda=%.2f$' % (i+1, lf)
  ax.plot(xm, am, c, ls='-', lw=2.0, label=lbl)
  lbl = r'$n=%i,\ \widetilde{\lambda}=%.2f$' % (i+1, lm)
  ax.plot(xf, af, c, ls='--',lw=2.0, label=lbl)

leg = ax.legend(loc='lower left', handlelength=3, fontsize='x-small')
leg.get_frame().set_alpha(0.0)
ax.set_ylim([-1.1,1.1])
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$u$')
ax.grid()
tight_layout()
savefig('../../images/fenics_intro/eigenvectors.pdf')

cn = array(cn, dtype='d')
af = zeros(1000)
am = zeros(n+1)
for i in range(n-1):
  af += cn[i]*Af[i]
  am += cn[i]*Am[i]

fig = figure(figsize=(5,3.5))
ax  = fig.add_subplot(111)

ax.plot(xf, af, 'k--', lw=2.0, label=r'$u(x)$ - Fourier')
ax.plot(xm, am, 'k-',  lw=2.0, label=r'$\widetilde{u}(x)$ - FEM')

leg = ax.legend(loc='center')
leg.get_frame().set_alpha(0.0)
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$u$')
ax.grid()
tight_layout()
savefig('../../images/fenics_intro/eigen_solution.pdf')



