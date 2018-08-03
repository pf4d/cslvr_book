from numpy import *
import matplotlib.pyplot as plt
import matplotlib        as mpl
import os
import re
import sys

mpl.rcParams['font.family']          = 'serif'
mpl.rcParams['legend.fontsize']      = 'small'
mpl.rcParams['text.usetex']          = True
mpl.rcParams['text.latex.preamble']  = ['\usepackage[mathscr]{euscript}']

# set the relavent directories :
#reg     = 'TV'
reg     = 'Tikhonov'
in_dir  = 'L_curve_results_' + reg + '/'
plt_dir = '../../images/data_assimilation/l_curve/'

# get the list of directories to plot :
alphas = []
Ds     = []
Js     = []
J1s    = []
J2s    = []
Rs     = []
ns     = []

for d in next(os.walk(in_dir))[1]:
  m = re.search('(alpha_)(\d\W\dE\W\d+)', d)
  if m is not None:
    do = in_dir + d + '/objective_ftnls_history/'
    alphas.append(float(m.group(2)))
    Ds.append(loadtxt(do  + 'Ds.txt'))
    Js.append(loadtxt(do  + 'Js.txt'))
    J1s.append(loadtxt(do + 'J1s.txt'))
    J2s.append(loadtxt(do + 'J2s.txt'))
    Rs.append(loadtxt(do  + 'Rs.txt'))
    ns.append(len(Js[-1]))
alphas = array(alphas) 
Ds     = array(Ds) 
Js     = array(Js) 
J1s    = array(J1s) 
J2s    = array(J2s) 
Rs     = array(Rs) 
ns     = array(ns)

idx    = argsort(alphas)
alphas = alphas[idx]
Ds     = Ds[idx]
Js     = Js[idx]
J1s    = J1s[idx]
J2s    = J2s[idx]
Rs     = Rs[idx]
ns     = ns[idx]
xn     = cumsum(ns - 1)

# define the colors we need :
#===============================================================================

cmap   = plt.get_cmap('viridis')
colors = [ cmap(x) for x in linspace(0, 1, len(alphas)) ]

dkred  = '#820000'
pink   = '#ff8e8e'
ltpink = '#ffbaba'


# plot functionals :
#===============================================================================

fig = plt.figure(figsize=(6,2.5))
ax  = fig.add_subplot(111)

k    = 0
ints = [0]
for i,c in zip(range(len(alphas)), colors):
  xi = arange(k, k + ns[i])
  ints.append(xi.max())
  if i == 0:
    #ax.plot(xi, Js[i], '-',  c='k', lw=2.0,
    #        label = r'$\mathscr{I}$')
    ax.plot(xi, J1s[i], '-',  c='0.5', lw=2.0,
            label = r'$\mathscr{I}_1$')
    ax.plot(xi, J2s[i], '-',  c='k', lw=2.0,
            label = r'$\mathscr{I}_2$')
    ax.plot(xi, Rs[i], '-', c='r', lw=2.0,
            label = r'$\mathscr{R}$')
  else:
    #ax.plot(xi, Js[i], '-',  c='k', lw=2.0)
    ax.plot(xi, J1s[i], '-',  c='0.5', lw=2.0)
    ax.plot(xi, J2s[i], '-',  c='k', lw=2.0)
    ax.plot(xi, Rs[i], '-', c='r', lw=2.0)
  k += ns[i] - 1
ints = array(ints)

ax.grid()
  
if reg == 'Tikhonov':
  gamma = '\gamma_3'
elif reg == 'TV':
  gamma = '\gamma_4'

label = []
for i in alphas:
  label.append(r'$%s = %g$' % (gamma, i))

ax.set_xticks(ints)
ax.set_xticklabels(label, size='small', ha='left')#, rotation=-45)
ax.set_xlabel(r'relative iteration')

ax.set_yscale('log')

leg = ax.legend(loc='upper left', ncol=3)
leg.get_frame().set_alpha(0.0)

plt.tight_layout()
plt.savefig(plt_dir + reg + '_convergence.pdf')
plt.close(fig)

# plot L-curve :
#===============================================================================

fin_Js = Js[:,-1]
fin_Rs = Rs[:,-1]

fig = plt.figure(figsize=(6,2.5))
ax  = fig.add_subplot(111)

ax.plot(fin_Js, fin_Rs, 'k-', lw=2.0)

ax.grid()

for i,c in zip(range(len(alphas)), colors):
  ax.plot(fin_Js[i], fin_Rs[i], 'o',  c=c, lw=2.0,
          label = r'$%s = %g$' % (gamma, alphas[i]))

Js_label = []
for i in fin_Js:
  Js_label.append(r'$%.1e$' % i)

ax.set_xlabel(r'$\mathscr{I}^*$')
ax.set_ylabel(r'$\mathscr{R}^*$')

leg = ax.legend(loc='upper right', ncol=2)
leg.get_frame().set_alpha(0.0)

ax.set_yscale('log')
#ax.set_xscale('log')

plt.tight_layout()
plt.savefig(plt_dir + reg + '_l_curve.pdf')
plt.close(fig)



