from numpy import *
import matplotlib.pyplot as plt
import matplotlib        as mpl
import os

mpl.rcParams['font.family']          = 'serif'
mpl.rcParams['legend.fontsize']      = 'xx-small'
mpl.rcParams['text.usetex']          = True
mpl.rcParams['text.latex.preamble']  = ['\usepackage[mathscr]{euscript}']

# set the relavent directories :
Fb_dir  = 'ps_results/Fb/'
ze_dir  = 'ps_results/zero_energy/'
plt_dir = '../../images/tmc/plane_strain/'

# plot convergence :
Fb_abs_err = loadtxt(Fb_dir + 'tmc/convergence_history/abs_err.txt')
Fb_tht_nrm = loadtxt(Fb_dir + 'tmc/convergence_history/theta_norm.txt')
ze_abs_err = loadtxt(ze_dir + 'tmc/convergence_history/abs_err.txt')
ze_tht_nrm = loadtxt(ze_dir + 'tmc/convergence_history/theta_norm.txt')

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


out = get_data(Fb_dir)

Is, Ds, Js, ns, xn = out

fig1 = plt.figure(figsize=(6,2.5))
ax1  = fig1.add_subplot(111)
ax2  = ax1.twinx()
k    = 0
ints = [0]
for i in range(len(Is)):
  xi = arange(k, k + ns[i])
  ints.append(xi.max())
  
  ax1.plot(xi, Js[i], '-', c='k', lw=2.0)
  ax2.plot(xi, Ds[i], '-', c='r', lw=2.0)

  k += ns[i] - 1
ints = array(ints)

ax1.grid()

ax1.set_xlim([0, max(xn.max(), xn.max())])
ax2.set_xlim([0, max(xn.max(), xn.max())])

ax1.set_xlabel('iteration')
ax1.set_ylabel(r'$\mathscr{J}$')
ax2.set_ylabel(r'$\Vert \theta - \theta_c \Vert_{\infty}$')

ax1.set_xticks(ints)
ax1.set_xticklabels(arange(len(ints)+1))
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useOffset=False)
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useOffset=False)
ax2.tick_params(axis='y', colors='r')

yloc1 = plt.MaxNLocator(4)
yloc2 = plt.MaxNLocator(4)
ax1.yaxis.set_major_locator(yloc1)

ax1.set_yscale('log')
ax2.set_yscale('log')
#ax3.set_yscale('log')
ax2.yaxis.get_offset_text().set_color('r')
ax2.yaxis.label.set_color('r')

plt.tight_layout()
plt.savefig(plt_dir + 'opt_convergence.pdf')
plt.close(fig1)


fig2 = plt.figure(figsize=(6,2.5))
ax3  = fig2.add_subplot(111)
ax4  = ax3.twinx()

dkred  = '#820000'
pink   = '#ff8e8e'
ltpink = '#ffbaba'

x_fb = arange(1, len(Fb_abs_err)+1)
x_ze = arange(1, len(ze_abs_err)+1)

ax3.plot(x_fb, Fb_abs_err, 'o-', c='k',   lw=2.0,
         label=r'$F_b^*$')
ax4.plot(x_fb, Fb_tht_nrm, 'o-', c='r',   lw=2.0,
         label=r'$F_b^*$')
ax3.plot(x_ze, ze_abs_err, 'o-', c='0.5', lw=2.0,
         label=r'$F_b = M_b \rho / \rho_w$')
ax4.plot(x_ze, ze_tht_nrm, 'o-', c=pink,   lw=2.0,
         label=r'$F_b = M_b \rho / \rho_w$')

leg1 = ax3.legend(loc='upper right')
leg2 = ax4.legend(loc='center right')
leg1.get_frame().set_alpha(0.0)
leg2.get_frame().set_alpha(0.0)

ax3.grid()
ax3.set_xlabel('iteration')
#ax4.set_xlabel('iteration')
ax3.set_ylabel(r'$\Vert \theta_n - \theta_{n-1} \Vert_2$')
ax4.set_ylabel(r'$\Vert \theta_n \Vert_2$')
ax3.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useOffset=False)
ax4.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useOffset=False)
ax4.tick_params(axis='y', colors='r')
yloc3 = plt.MaxNLocator(4)
yloc4 = plt.MaxNLocator(4)
ax3.yaxis.set_major_locator(yloc3)
ax4.yaxis.set_major_locator(yloc4)
ax4.yaxis.get_offset_text().set_color('r')
ax4.yaxis.label.set_color('r')

plt.tight_layout()
plt.savefig(plt_dir + 'tmc_convergence.pdf')
plt.close(fig2)



