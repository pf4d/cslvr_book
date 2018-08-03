from pylab import *

mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['legend.fontsize'] = 'medium'

rhoi = 917.0
rhow = 1000.0
ki   = 2.1
kw   = 0.561
ci   = 2009.0
cw   = 4217.6
a    = 60*60*24*365
     
h1   = 1000.0
h2   = 100000.0
u    = logspace(log(5e-8), log(100000.0),100000, base=e)
#u    = linspace(0, 100000.0, 1000)

ks   = ki
cs   = ci
rhos = rhoi

kap  = a * ks/(cs*rhos)

A1  = h1*u / (2*kap)
xi1 = 1/tanh(A1) - 1/A1

A2  = h2*u / (2*kap)
xi2 = 1/tanh(A2) - 1/A2

fig = figure(figsize=(8,7))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

ax1.plot(u, A2, 'r-', lw=2.0, label=r'$h = 100$ km')
ax1.plot(u, A1, 'k-', lw=2.0, label=r'$h = 1$ km')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim([u.min(),  u.max()])
ax1.set_ylim([A1.min(), A2.max()])
ax1.set_xlabel(r'$\Vert \mathbf{u} \Vert$')
ax1.set_ylabel(r'$P_e$')
ax1.legend(loc='lower right')
ax1.grid()

ax2.plot(u, xi2, 'r-', lw=2.0, label=r'$h = 100$ km')
ax2.plot(u, xi1, 'k-', lw=2.0, label=r'$h = 1$ km')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim([u.min(), u.max()])
ax2.set_ylim([xi1.min(), 10])
ax2.set_xlabel(r'$\Vert \mathbf{u} \Vert$')
ax2.set_ylabel(r'$\xi$')
ax2.legend(loc='lower right')
ax2.grid()

tight_layout()
savefig('../../images/internal_energy/alpha_new.pdf')
show()



