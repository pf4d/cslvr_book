from pylab import *
from matplotlib import colors
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
  '''
  Function to offset the "center" of a colormap. Useful for
  data with a negative min and positive max and you want the
  middle of the colormap's dynamic range to be at zero

  Input
  -----
    cmap : The matplotlib colormap to be altered
    start : Offset from lowest point in the colormap's range.
        Defaults to 0.0 (no lower ofset). Should be between
        0.0 and `midpoint`.
    midpoint : The new center of the colormap. Defaults to 
        0.5 (no shift). Should be between 0.0 and 1.0. In
        general, this should be  1 - vmax/(vmax + abs(vmin))
        For example if your data range from -15.0 to +5.0 and
        you want the center of the colormap at 0.0, `midpoint`
        should be set to  1 - 5/(5 + 15)) or 0.75
    stop : Offset from highets point in the colormap's range.
        Defaults to 1.0 (no upper ofset). Should be between
        `midpoint` and 1.0.
  '''
  cdict = {
      'red': [],
      'green': [],
      'blue': [],
      'alpha': []
  }

  # regular index to compute the colors
  reg_index = np.linspace(start, stop, 257)

  # shifted index to match the data
  shift_index = np.hstack([
    np.linspace(0.0, midpoint, 128, endpoint=False), 
    np.linspace(midpoint, 1.0, 129, endpoint=True)
  ])

  for ri, si in zip(reg_index, shift_index):
    r, g, b, a = cmap(ri)

    cdict['red'].append((si, r, r))
    cdict['green'].append((si, g, g))
    cdict['blue'].append((si, b, b))
    cdict['alpha'].append((si, a, a))

  newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
  plt.register_cmap(cmap=newcmap)

  return newcmap

mpl.rcParams['font.family']          = 'serif'
mpl.rcParams['legend.fontsize']      = 'medium'
mpl.rcParams['text.usetex']          = True
mpl.rcParams['text.latex.preamble']  = ['\usepackage[mathscr]{euscript}']

rhoi  = 910.0
rhow  = 1000.0
g     = 9.81
gamma = 9.8e-8 
Tw    = 273.15
ki    = 2.1
kw    = 0.561
ci    = 2009.0
cw    = 4217.6
a     = 31556926.0
wmax  = 0.01
n     = 3.0
R     = 8.3144621
k     = 1000

H     = linspace(0, 5000, k)
w     = linspace(0, wmax, k)
t     = linspace(-50, 0,  k)
tw    = linspace(-8, 0,  k)
a_T   = zeros(k)
Q_T   = zeros(k)
a_T_l = 3.985e-13 * a
a_T_u = 1.916e3 * a
Q_T_l = 6e4
Q_T_u = 13.9e4
Tc    = - 10

a_T[t <= Tc] = a_T_l
a_T[t >  Tc] = a_T_u
Q_T[t <= Tc] = Q_T_l
Q_T[t >  Tc] = Q_T_u

W,T = meshgrid(w,tw)

a_T_m = zeros(shape(W))
Q_T_m = zeros(shape(W))

a_T_m[T <= Tc] = a_T_l
a_T_m[T >  Tc] = a_T_u
Q_T_m[T <= Tc] = Q_T_l
Q_T_m[T >  Tc] = Q_T_u

A  = a_T*exp(-Q_T/(R*(t + Tw)))
Aw = a_T_m*(1 + 181.25*W)*exp(-Q_T_m/(R*(T+Tw)))

fig = figure(figsize=(7,4))
ax  = gridspec.GridSpec(1,2, width_ratios=[1,1.25])
ax1 = plt.subplot(ax[0,0])
ax2 = plt.subplot(ax[0,1])

spTick = a_T_l*exp(-Q_T_l/(R*(Tc + Tw)))
ax1.set_yscale('log', basey=e)
ax1.set_yticks([A.min(), 1e-18, spTick, A.max()])

ax1.plot(t, A, 'k', lw=2.0, label=r'$A$')
ax1.set_xlim([t.min(), t.max()])
ax1.set_ylim([A.min(), A.max()])
ax1.set_xticks([-50, -40, -30, -20, -10, 0])
ax1.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%i'))
ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1e'))
ax1.set_xlabel(r"$T'$")
ax1.set_ylabel(r'$A$')
ax1.grid()

lt = [(-7, 1.5e-3),
      (-5.2, 4e-3),
      (-4.5, 5.3e-3),
      (-3.2, 6.5e-3),
      (-2.6, 7.2e-3),
      (-2, 7.8e-3),
      (-1.5, 8.2e-3),
      (-1.4, 8.5e-3),
      (-0.8, 8.8e-3)]

cm  = get_cmap('gist_yarg')
#lv1 = array([4e-17, 6e-17, 8e-17, 1e-16, 1.1e-16, 1.2e-16, 3e-16])
lv1 = linspace(4e-17, 4e-16, 10)
cs2 = ax2.contour(T,W,Aw, linewidths=2.0, colors='k', levels=lv1)
ax2.clabel(cs2, inline=1, colors='k', fmt='%.1e', manual=lt)
ax2.set_xlabel(r"$T'$")
ax2.set_ylabel(r'$W$')
ax2.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0e'))
ax2.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%i'))
ax2.grid()

tight_layout()
savefig('../images/rate_factor_new.pdf')
close(fig)




