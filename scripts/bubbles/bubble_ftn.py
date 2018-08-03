from pylab import *

x = linspace(0,1,1000)

phi_1 = 1 - x
phi_2 = x

mpl.rcParams['font.family']          = 'serif'
mpl.rcParams['text.usetex']          = True
mpl.rcParams['text.latex.preamble']  = ['\usepackage{fouriernc}']

fig = figure(figsize=(5,3.5))
ax  = fig.add_subplot(111)

ax.set_ylim(0.0, 1.10)

ax.plot(x, phi_1,         'k-',  lw=2.0, label=r"$\psi_1$")
ax.plot(x, phi_2,         'k--', lw=2.0, label=r"$\psi_2$")
ax.plot(x, 4*phi_1*phi_2, 'r',   lw=2.0, label=r"$\phi'$")

leg = ax.legend(loc='lower center')
leg.get_frame().set_alpha(0.0)
ax.set_xlabel(r"$x$")
ax.grid()
tight_layout()
savefig("../../images/bubbles/bubble_new.pdf")
