from pylab import *

x = linspace(0,1,1000)

u1 = 270/31.
u2 = -140/31.

u = 1 + u1*(x**2 - 2*x) + u2*(x**3 - 3*x)
e = cosh(x) - tanh(1)*sinh(x)

mpl.rcParams['font.family']          = 'serif'
mpl.rcParams['text.usetex']          = True
mpl.rcParams['text.latex.preamble']  = ['\usepackage{fouriernc}']

fig = figure(figsize=(5,3))
ax  = fig.add_subplot(111)

ax.plot(x, u, 'k',   lw=2.0, label='approximate')
ax.plot(x, e, 'k--', lw=2.0, label='exact')

leg = ax.legend(loc="upper left")
leg.get_frame().set_alpha(0.0)
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$u$')
ax.grid()

tight_layout()
savefig("../../images/fenics_intro/weight_int.pdf")


