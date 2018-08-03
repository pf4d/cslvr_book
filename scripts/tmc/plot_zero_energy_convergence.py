from numpy import *
import matplotlib.pyplot as plt
import matplotlib        as mpl
import os

# these control the basal-boundary condition :
e_mode  = 'zero_energy'

mpl.rcParams['font.family']          = 'serif'
mpl.rcParams['legend.fontsize']      = 'xx-small'
mpl.rcParams['text.usetex']          = True
mpl.rcParams['text.latex.preamble']  = ['\usepackage[mathscr]{euscript}']

# set the relavent directories :
in_dir  = 'ps_results/' + e_mode + '/'
plt_dir = '../../images/tmc/plane_strain/' + e_mode + '/'

# plot convergence :
abs_err = loadtxt(in_dir + 'tmc/convergence_history/abs_err.txt')
tht_nrm = loadtxt(in_dir + 'tmc/convergence_history/theta_norm.txt')

fig = plt.figure(figsize=(6,2.3))
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()

x = arange(1,len(abs_err)+1)

ax1.plot(x, abs_err, 'o-', c='k', lw=2.0)
ax2.plot(x, tht_nrm, 'o-', c='r', lw=2.0)

ax1.grid()

#ax1.set_xlim([0, max(xn.max(), xn.max())])
#ax2.set_xlim([0, max(xn.max(), xn.max())])

ax1.set_xlabel('iteration')
ax1.set_ylabel(r'$\Vert \theta_n - \theta_{n-1} \Vert_2$')
ax2.set_ylabel(r'$\Vert \theta_n \Vert_2$')

ax1.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useOffset=False)
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useOffset=False)
ax2.tick_params(axis='y', colors='r')

yloc1 = plt.MaxNLocator(4)
ax1.yaxis.set_major_locator(yloc1)

ax2.yaxis.get_offset_text().set_color('r')
ax2.yaxis.label.set_color('r')

plt.tight_layout()
plt.savefig(plt_dir + 'convergence.pdf')
plt.close(fig)




