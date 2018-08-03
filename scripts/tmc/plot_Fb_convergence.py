from numpy import *
import matplotlib.pyplot as plt
import matplotlib        as mpl
import os

e_mode  = 'Fb'

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

def get_data(direc):
  Ds = []
  Js = []
  Is = []
  ns = []
  
  for d in next(os.walk(direc + 'tmc/'))[1]:
    try:
      i = int(d)
      dn = direc + 'tmc/' + d + '/objective_ftnls_history/'
      Is.append(i)
      Ds.append(loadtxt(dn + 'Ds.txt'))
      Js.append(loadtxt(dn + 'Js.txt'))
      ns.append(len(Js[-1]))
    except ValueError:
      pass
  
  Ds = array(Ds)
  Js = array(Js)
  Is = array(Is)
  ns = array(ns)
  
  idx = argsort(Is)
  Ds  = Ds[idx]
  Js  = Js[idx]
  ns  = ns[idx]
  xn  = cumsum(ns - 1)

  return (Is, Ds, Js, ns, xn)


out = get_data(in_dir)

Is, Ds, Js, ns, xn = out

Jmax = Js.max()
amax = abs_err.max()
tmax = tht_nrm.max()

fig = plt.figure(figsize=(6,4))
ax1 = fig.add_subplot(211)
ax2 = ax1.twinx()
ax3 = fig.add_subplot(212)
ax4 = ax3.twinx()

#ax1 = fig.add_subplot(411)
#ax3 = fig.add_subplot(413)
#ax4 = fig.add_subplot(414)

k = 0

ax3.plot(xn, abs_err, 'o-', c='k', lw=2.0)
ax4.plot(xn, tht_nrm, 'o-', c='r', lw=2.0)


for i in range(len(Is)):
  xi = arange(k, k + ns[i])
  
  ax1.plot(xi, Js[i], '-', c='k', lw=2.0)
  ax2.plot(xi, Ds[i], '-', c='r', lw=2.0)

  k += ns[i] - 1

ax1.grid()
ax3.grid()

ax1.set_xlim([0, max(xn.max(), xn.max())])
ax2.set_xlim([0, max(xn.max(), xn.max())])
ax3.set_xlim([0, max(xn.max(), xn.max())])
ax4.set_xlim([0, max(xn.max(), xn.max())])

ax3.set_ylim([0.0, abs_err.max()])

ax1.set_xticklabels([])

#ax1.set_xlabel('iteration')
ax3.set_xlabel('iteration')
#ax4.set_xlabel('iteration')
ax1.set_ylabel(r'$\mathscr{J}$')
ax2.set_ylabel(r'$\Vert \theta - \theta_c \Vert_{\infty}$')
ax3.set_ylabel(r'$\Vert \theta_n - \theta_{n-1} \Vert_2$')
ax4.set_ylabel(r'$\Vert \theta_n \Vert_2$')

ax1.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useOffset=False)
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useOffset=False)
ax3.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useOffset=False)
ax4.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useOffset=False)
ax4.tick_params(axis='y', colors='r')
ax2.tick_params(axis='y', colors='r')

yloc1 = plt.MaxNLocator(4)
yloc2 = plt.MaxNLocator(4)
yloc3 = plt.MaxNLocator(4)
yloc4 = plt.MaxNLocator(4)
ax1.yaxis.set_major_locator(yloc1)
ax3.yaxis.set_major_locator(yloc3)
ax4.yaxis.set_major_locator(yloc4)

ax1.set_yscale('log')
ax2.set_yscale('log')
#ax3.set_yscale('log')
ax2.yaxis.get_offset_text().set_color('r')
ax2.yaxis.label.set_color('r')
ax4.yaxis.get_offset_text().set_color('r')
ax4.yaxis.label.set_color('r')

plt.tight_layout()
plt.savefig(plt_dir + 'convergence.pdf')
plt.close(fig)




