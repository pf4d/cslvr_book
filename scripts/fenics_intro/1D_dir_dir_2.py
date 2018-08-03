from pylab  import *
from fenics import *

mesh = IntervalMesh(1000,-1,1)
Q    = FunctionSpace(mesh, 'CG', 1)

u    = interpolate(Expression('pow(x[0], 6)'),   Q)
v    = interpolate(Expression('pow(x[0], 2)'),   Q)
dudv = interpolate(Expression('3*pow(x[0], 4)'), Q)

dudv_1 = u.dx(0) * 1/v.dx(0)

x    = mesh.coordinates()[:,0][::-1]
u_v  = u.vector().array()
v_v  = v.vector().array()
d_va = dudv.vector().array()
d_v1 = project(dudv_1).vector().array()

mpl.rcParams['font.family']          = 'serif'
mpl.rcParams['text.usetex']          = True
mpl.rcParams['text.latex.preamble']  = ['\usepackage{fouriernc}']

fig  = figure(figsize=(5,3.5))
ax   = fig.add_subplot(111)

purp = '#880cbc'
grun = '#77f343'

ax.plot(x, u_v,  'k',   lw=2.0, label=r'$u(x)=x^6$')
ax.plot(x, v_v,  'r',   lw=2.0, label=r'$v(x)=x^2$')

ax.plot(x, d_va, color=grun, ls='-',  lw=2.0, 
        label=r'$\frac{du}{dv} = 3x^4$')
ax.plot(x, d_v1, color=purp, ls='--', lw=2.0,
        label=r'$\frac{du}{dv}$ - FEniCS')

ax.grid()
ax.set_xlabel(r'$x$')
leg = ax.legend(loc='upper center', ncol=2, fontsize='medium')
leg.get_frame().set_alpha(0.0)

tight_layout()
savefig('../../images/fenics_intro/1D_dir_dir_2.pdf')



