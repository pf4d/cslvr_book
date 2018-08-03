from cslvr    import *
from fenics   import *
import matplotlib.pyplot as plt
import numpy             as np
import sys

# set the relavent directories :
base_dir = 'dump/jakob_small/inversion_Wc_0.01/10/'
in_dir   = base_dir
out_dir  = '../../images/internal_energy/jakob_results/inversion_Wc_0.01/'
var_dir  = 'dump/vars_jakobshavn_small/'

if not os.path.exists(out_dir):
  os.makedirs(out_dir)

# laod the HDF5 files :
fdata   = HDF5File(mpi_comm_world(), var_dir + 'state.h5',       'r')
f       = HDF5File(mpi_comm_world(), in_dir  + 'inverted_10.h5', 'r')

# not deformed mesh :
mesh    = Mesh('dump/meshes/jakobshavn_3D_small_block.xml.gz')

# create 3D model for stokes solves :
model = D3Model(mesh=mesh, out_dir=out_dir)

# initialize the variables :
model.init_S(fdata)
model.init_B(fdata)
model.init_p(f)
model.init_theta(f)
model.init_T(f)
model.init_W(f)

#===============================================================================
drg  = DataFactory.get_rignot()

bc = '#880cbc'

lat_1 = 69.210
lat_2 = 69.168
lon_1 = -48.78
lon_2 = -48.759

dlat  = (lat_2 - lat_1) / 2.0
dlon  = (lon_2 - lon_1) / 2.0

lat_3 = lat_1 + dlat
lon_3 = lon_1 + dlon

lat_a   = [ 69.235,    lat_1, lat_2, lat_3]
lon_a   = [-48.686944, lon_1, lon_2, lon_3]
color_a = ['c', 'y', 'g',  bc]
style_a = ['o', 'p', '^', 's']

gamma = model.gamma(0)
Tw    = model.T_w(0)
L     = model.L(0)
a     = 146.3
b     = 7.253

#===============================================================================
# plot profiles :

x_a, y_a  = drg['pyproj_Proj'](lon_a, lat_a)

zmin = mesh.coordinates()[:,2].min()
zmax = mesh.coordinates()[:,2].max()

z_s = linspace(zmin, zmax, 100)

T_a = []
W_a = []
z_a = []

for x_w, y_w in zip(x_a, y_a):
  S    = model.S(x_w, y_w, 1.0)
  B    = model.B(x_w, y_w, 1.0)
  T_z  = []
  W_z  = []
  for z_w in z_s:
    theta_i = model.theta(x_w, y_w, z_w)
    p_i     = model.p(x_w, y_w, z_w)
    Tm_i    = Tw - gamma*p_i
    theta_m = a*Tm_i + b/2*Tm_i**2
    if theta_i > theta_m:
      W_z.append( (theta_i - theta_m)/L )
      T_z.append( Tm_i )
    else:
      W_z.append( 0.0 )
      T_z.append( (-a + np.sqrt(a**2 + 2*b*theta_i)) / b )
  
  T_z = array(T_z)
  W_z = array(W_z)

  # get z-coordinates :  
  z_z = []
  for z_w in z_s:
    z_i = (z_w / zmax) - 1# * (S - B) - (S - B)
    z_z.append(z_i)
  z_z = array(z_z)
  
  T_a.append(T_z)
  W_a.append(W_z)
  z_a.append(z_z)

T_a = array(T_a)
W_a = array(W_a)
z_a = array(z_a)

if not os.path.exists(base_dir + 'profile_data'):
  os.makedirs(base_dir + 'profile_data')

np.savetxt(base_dir + 'profile_data/T.txt', T_a)
np.savetxt(base_dir + 'profile_data/W.txt', W_a)
np.savetxt(base_dir + 'profile_data/z.txt', z_a)

fig = figure(figsize=(5,5))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

for T_i, W_i, z_i, color_i, style_i in zip(T_a, W_a, z_a, color_a, style_a):
  ax1.plot(T_i, z_i, color=color_i, lw=3.0)
  ax2.plot(W_i, z_i, color=color_i, lw=3.0)

# get every 10th value to plot nodes of mesh :
for T_i, W_i, z_i, color_i, style_i in zip(T_a, W_a, z_a, color_a, style_a):
  ax1.plot(T_i[::10], z_i[::10], marker=style_i, color=color_i, ls='None')
  ax2.plot(W_i[::10], z_i[::10], marker=style_i, color=color_i, ls='None')


ax2.set_yticklabels([])
#ax2.set_xlim([0,0.15])

xloc1 = plt.MaxNLocator(4)
xloc2 = plt.MaxNLocator(4)
ax1.xaxis.set_major_locator(xloc1)
ax2.xaxis.set_major_locator(xloc2)

ax2.ticklabel_format(axis='x', style='sci', scilimits=(0,0), useOffset=False)

ax1.set_xlabel(r'$T$')
ax2.set_xlabel(r'$W$')
ax1.set_ylabel(r'relative depth')
#ax2.tick_params(axis='x', colors='r')
#ax2.xaxis.label.set_color('r')
ax1.grid()
ax2.grid()
plt.tight_layout()
plt.savefig(out_dir + 'profile_plot.pdf')
plt.close(fig)



