from cslvr    import *
from fenics   import *
import matplotlib.pyplot as plt
import numpy             as np
import sys

# set the relavent directories :
base_dir = 'dump/jakob_small/inversion_Wc_0.01/10/'
in_dir   = base_dir
out_dir  = '../../images/internal_energy/jakob_results/inversion_Wc_0.01/stress_balance/'

# not deformed mesh :
mesh    = Mesh('dump/meshes/jakobshavn_3D_small_block.xml.gz')

# create 3D model for stokes solves :
d3model = D3Model(mesh, out_dir)

#===============================================================================
# retrieve the bed mesh :
d3model.form_bed_mesh()
d3model.form_srf_mesh()

# create 2D model for balance velocity :
bedmodel = D2Model(d3model.bedmesh, out_dir)
srfmodel = D2Model(d3model.srfmesh, out_dir)

#===============================================================================
# open the hdf5 file :
f     = HDF5File(mpi_comm_world(), in_dir  + 'stress.h5', 'r')

# initialize the variables :
d3model.init_N_ii(f)
d3model.init_N_ij(f)
d3model.init_N_iz(f)
d3model.init_N_ji(f)
d3model.init_N_jj(f)
d3model.init_N_jz(f)
d3model.init_N_zi(f)
d3model.init_N_zj(f)
d3model.init_N_zz(f)
d3model.init_M_ii(f)
d3model.init_M_ij(f)
d3model.init_M_iz(f)
d3model.init_M_ji(f)
d3model.init_M_jj(f)
d3model.init_M_jz(f)
d3model.init_M_zi(f)
d3model.init_M_zj(f)
d3model.init_M_zz(f)

# 2D model gets balance-velocity appropriate variables initialized :
srfmodel.assign_submesh_variable(srfmodel.M_ii, d3model.M_ii)
srfmodel.assign_submesh_variable(srfmodel.M_ij, d3model.M_ij)
srfmodel.assign_submesh_variable(srfmodel.M_iz, d3model.M_iz)
srfmodel.assign_submesh_variable(srfmodel.M_ji, d3model.M_ji)
srfmodel.assign_submesh_variable(srfmodel.M_jj, d3model.M_jj)
srfmodel.assign_submesh_variable(srfmodel.M_jz, d3model.M_jz)
srfmodel.assign_submesh_variable(srfmodel.M_zi, d3model.M_zi)
srfmodel.assign_submesh_variable(srfmodel.M_zj, d3model.M_zj)
srfmodel.assign_submesh_variable(srfmodel.M_zz, d3model.M_zz)
srfmodel.assign_submesh_variable(srfmodel.N_ii, d3model.N_ii)
srfmodel.assign_submesh_variable(srfmodel.N_ij, d3model.N_ij)
srfmodel.assign_submesh_variable(srfmodel.N_iz, d3model.N_iz)
srfmodel.assign_submesh_variable(srfmodel.N_ji, d3model.N_ji)
srfmodel.assign_submesh_variable(srfmodel.N_jj, d3model.N_jj)
srfmodel.assign_submesh_variable(srfmodel.N_jz, d3model.N_jz)
srfmodel.assign_submesh_variable(srfmodel.N_zi, d3model.N_zi)
srfmodel.assign_submesh_variable(srfmodel.N_zj, d3model.N_zj)
srfmodel.assign_submesh_variable(srfmodel.N_zz, d3model.N_zz)

#===============================================================================
drg  = DataFactory.get_rignot()

bc      = '#880cbc'

#lat_1 = 69.210
#lat_2 = 69.168
#lon_1 = -48.78
#lon_2 = -48.759
#
#dlat  = (lat_2 - lat_1) / 2.0
#dlon  = (lon_2 - lon_1) / 2.0
#
#lat_3 = lat_1 + dlat
#lon_3 = lon_1 + dlon
#
#lat_a   = [ 69.235,    lat_1, lat_2, lat_3]
#lon_a   = [-48.686944, lon_1, lon_2, lon_3]
#color_a = ['c', 'y', 'g',  bc]
#style_a = ['o', 'p', '^', 's']
#
#plot_pts = {'lat'   : lat_a,
#            'lon'   : lon_a,
#            'style' : style_a,
#            'color' : color_a}
#
#zoom_box_kwargs = {'zoom'             : 5.8,      # ammount to zoom 
#                   'loc'              : 1,        # location of box
#                   'loc1'             : 2,        # loc of first line
#                   'loc2'             : 3,        # loc of second line
#                   'x1'               : 40000,    # first x-coord
#                   'y1'               : 80000,    # first y-coord
#                   'x2'               : 90000,    # second x-coord
#                   'y2'               : 105000,   # second y-coord
#                   'scale_font_color' : bc,       # scale font color
#                   'scale_length'     : 20,       # scale length in km
#                   'scale_loc'        : 1,        # 1=top, 2=bottom
#                   'plot_grid'        : True,     # plot the triangles
#                   'axes_color'       : bc,       # color of axes
#                   'plot_points'      : plot_pts} # dict of points to plot

params = {'llcrnrlat'    : 68.99,
          'urcrnrlat'    : 69.31,
          'llcrnrlon'    : -49.8,
          'urcrnrlon'    : -48.3,
          'scale_color'  : bc,
          'scale_length' : 50,
          'scale_loc'    : 1,
          'figsize'      : (7,4),
          'lat_interval' : 0.05,
          'lon_interval' : 0.25,
          'plot_grid'    : False,
          'plot_scale'   : False,
          'axes_color'   : 'r'}


M_ii_min = srfmodel.M_ii.vector().min()
M_ii_max = srfmodel.M_ii.vector().max()

M_ij_min = srfmodel.M_ij.vector().min()
M_ij_max = srfmodel.M_ij.vector().max()

M_iz_min = srfmodel.M_iz.vector().min()
M_iz_max = srfmodel.M_iz.vector().max()

M_ji_min = srfmodel.M_ji.vector().min()
M_ji_max = srfmodel.M_ji.vector().max()

M_jj_min = srfmodel.M_jj.vector().min()
M_jj_max = srfmodel.M_jj.vector().max()

M_jz_min = srfmodel.M_jz.vector().min()
M_jz_max = srfmodel.M_jz.vector().max()

M_zi_min = srfmodel.M_zi.vector().min()
M_zi_max = srfmodel.M_zi.vector().max()

M_zj_min = srfmodel.M_zj.vector().min()
M_zj_max = srfmodel.M_zj.vector().max()

M_zz_min = srfmodel.M_zz.vector().min()
M_zz_max = srfmodel.M_zz.vector().max()


N_ii_min = srfmodel.N_ii.vector().min()
N_ii_max = srfmodel.N_ii.vector().max()

N_ij_min = srfmodel.N_ij.vector().min()
N_ij_max = srfmodel.N_ij.vector().max()

N_iz_min = srfmodel.N_iz.vector().min()
N_iz_max = srfmodel.N_iz.vector().max()

N_ji_min = srfmodel.N_ji.vector().min()
N_ji_max = srfmodel.N_ji.vector().max()

N_jj_min = srfmodel.N_jj.vector().min()
N_jj_max = srfmodel.N_jj.vector().max()

N_jz_min = srfmodel.N_jz.vector().min()
N_jz_max = srfmodel.N_jz.vector().max()

N_zi_min = srfmodel.N_zi.vector().min()
N_zi_max = srfmodel.N_zi.vector().max()

N_zj_min = srfmodel.N_zj.vector().min()
N_zj_max = srfmodel.N_zj.vector().max()

N_zz_min = srfmodel.N_zz.vector().min()
N_zz_max = srfmodel.N_zz.vector().max()



M_ii_lvls = np.array([M_ii_min, -8e4, -4e4, -2e4, -5e3, 
                        5e3, 2e4, 4e4, 8e4, M_ii_max])
M_ij_lvls = np.array([M_ij_min, -1.5e5, -5e4, -2.5e4, -5e3, 
                        5e3, 2.5e4, 5e4, 1.5e5, M_ij_max])
M_iz_lvls = np.array([M_iz_min, -8e4, -4e4, -2e4, -5e3, 
                        5e3, 2e4, 4e4, 8e4, M_iz_max])

M_ji_lvls = np.array([M_ji_min, -8e4, -4e4, -2e4, -5e3, 
                        5e3, 2e4, 4e4, 8e4, M_ji_max])
M_jj_lvls = np.array([M_jj_min, -1.5e5, -5e4, -2.5e4, -5e3, 
                        5e3, 2.5e4, 5e4, 1.5e5, M_jj_max])
M_jz_lvls = np.array([M_jz_min, -2e5, -1e5, -5e4, -5e3, 
                        5e3, 5e4, 1e5, 2e5, M_jz_max])

M_zi_lvls = np.array([M_zi_min, -8e4, -4e4, -2e4, -5e3, 
                        5e3, 2e4, 4e4, 8e4, M_zi_max])
M_zj_lvls = np.array([M_zj_min, -1.5e5, -4e4, -2e4, -5e3, 
                        5e3, 2e4, 4e4, 8e4, M_zj_max])
M_zz_lvls = np.array([M_zz_min, -2e5, -1e5, -5e4, -5e3, 
                        5e3, 5e4, 1e5, 2e5, M_zz_max])



N_ii_lvls = np.array([N_ii_min, -2e8, -1e8, -5e7, -1e7, 
                        1e7, 5e7, 1e8, 2e8, N_ii_max])
N_ij_lvls = np.array([N_ij_min, -2e8, -1e8, -5e7, -1e7, 
                        1e7, 5e7, 1e8, 2e8, N_ij_max])
N_iz_lvls = np.array([N_iz_min, -2e8, -1e8, -5e7, -1e7, 
                        1e7, 5e7, 1e8, 2e8, N_iz_max])

N_ji_lvls = np.array([N_ji_min, -2e8, -1e8, -5e7, -1e7, 
                        1e7, 5e7, 1e8, 2e8, N_ji_max])
N_jj_lvls = np.array([N_jj_min, -1.5e8, -1e8, -5e7, -1e7, 
                        1e7, 5e7, 1e8, 2e8, N_jj_max])
N_jz_lvls = np.array([N_jz_min, -2e8, -1e8, -5e7, -1e7, 
                        1e7, 5e7, 1e8, 2e8, N_jz_max])

N_zi_lvls = np.array([N_zi_min, -2e8, -1e8, -5e7, -1e7, 
                        1e7, 5e7, 1e8, 2e8, N_zi_max])
N_zj_lvls = np.array([N_zj_min, -2e8, -1e8, -5e7, -1e7, 
                        1e7, 5e7, 1e8, 2e8, N_zj_max])
N_zz_lvls = np.array([N_zz_min, -2e8, -1e8, -5e7, -1e7, 
                        1e7, 5e7, 1e8, 2e8, N_zz_max])


#===============================================================================
# plot :

#cmap = 'RdGy'
#cmap = 'viridis'
#cmap = 'inferno'
#cmap = 'plasma'
#cmap = 'magma'
cmap = 'gist_yarg'
  

plotIce(drg, srfmodel.M_ii, name='M_ii', direc=out_dir,
        cmap='RdGy',  scale='lin',
        levels=M_ii_lvls, tp=False, tpAlpha=0.2,
        extend='neither', show=False, ext='.pdf',
        params=params)

plotIce(drg, srfmodel.M_ij, name='M_ij', direc=out_dir,
        cmap='RdGy',  scale='lin',
        levels=M_ij_lvls, tp=False, tpAlpha=0.2,
        extend='neither', show=False, ext='.pdf',
        params=params)

plotIce(drg, srfmodel.M_iz, name='M_iz', direc=out_dir,
        cmap='RdGy',  scale='lin',
        levels=M_iz_lvls, tp=False, tpAlpha=0.2,
        extend='neither', show=False, ext='.pdf',
        params=params)

  
plotIce(drg, srfmodel.M_ji, name='M_ji', direc=out_dir,
        cmap='RdGy',  scale='lin',
        levels=M_ji_lvls, tp=False, tpAlpha=0.2,
        extend='neither', show=False, ext='.pdf',
        params=params)

plotIce(drg, srfmodel.M_jj, name='M_jj', direc=out_dir,
        cmap='RdGy',  scale='lin',
        levels=M_jj_lvls, tp=False, tpAlpha=0.2,
        extend='neither', show=False, ext='.pdf',
        params=params)

plotIce(drg, srfmodel.M_jz, name='M_jz', direc=out_dir,
        cmap='RdGy',  scale='lin',
        levels=M_jz_lvls, tp=False, tpAlpha=0.2,
        extend='neither', show=False, ext='.pdf',
        params=params)

  
plotIce(drg, srfmodel.M_zi, name='M_zi', direc=out_dir,
        cmap='RdGy',  scale='lin',
        levels=M_zi_lvls, tp=False, tpAlpha=0.2,
        extend='neither', show=False, ext='.pdf',
        params=params)

plotIce(drg, srfmodel.M_zj, name='M_zj', direc=out_dir,
        cmap='RdGy',  scale='lin',
        levels=M_zj_lvls, tp=False, tpAlpha=0.2,
        extend='neither', show=False, ext='.pdf',
        params=params)

plotIce(drg, srfmodel.M_zz, name='M_zz', direc=out_dir,
        cmap='RdGy',  scale='lin',
        levels=M_zz_lvls, tp=False, tpAlpha=0.2,
        extend='neither', show=False, ext='.pdf',
        params=params)
 

  
plotIce(drg, srfmodel.N_ii, name='N_ii', direc=out_dir,
        cmap='RdGy',  scale='lin',
        levels=N_ii_lvls, tp=False, tpAlpha=0.2,
        extend='neither', show=False, ext='.pdf',
        params=params)

plotIce(drg, srfmodel.N_ij, name='N_ij', direc=out_dir,
        cmap='RdGy',  scale='lin',
        levels=N_ij_lvls, tp=False, tpAlpha=0.2,
        extend='neither', show=False, ext='.pdf',
        params=params)

plotIce(drg, srfmodel.N_iz, name='N_iz', direc=out_dir,
        cmap='RdGy',  scale='lin',
        levels=N_iz_lvls, tp=False, tpAlpha=0.2,
        extend='neither', show=False, ext='.pdf',
        params=params)

  
plotIce(drg, srfmodel.N_ji, name='N_ji', direc=out_dir,
        cmap='RdGy',  scale='lin',
        levels=N_ji_lvls, tp=False, tpAlpha=0.2,
        extend='neither', show=False, ext='.pdf',
        params=params)

plotIce(drg, srfmodel.N_jj, name='N_jj', direc=out_dir,
        cmap='RdGy',  scale='lin',
        levels=N_jj_lvls, tp=False, tpAlpha=0.2,
        extend='neither', show=False, ext='.pdf',
        params=params)

plotIce(drg, srfmodel.N_jz, name='N_jz', direc=out_dir,
        cmap='RdGy',  scale='lin',
        levels=N_jz_lvls, tp=False, tpAlpha=0.2,
        extend='neither', show=False, ext='.pdf',
        params=params)

  
plotIce(drg, srfmodel.N_zi, name='N_zi', direc=out_dir,
        cmap='RdGy',  scale='lin',
        levels=N_ii_lvls, tp=False, tpAlpha=0.2,
        extend='neither', show=False, ext='.pdf',
        params=params)

plotIce(drg, srfmodel.N_zj, name='N_zj', direc=out_dir,
        cmap='RdGy',  scale='lin',
        levels=N_ij_lvls, tp=False, tpAlpha=0.2,
        extend='neither', show=False, ext='.pdf',
        params=params)

plotIce(drg, srfmodel.N_zz, name='N_zz', direc=out_dir,
        cmap='RdGy',  scale='lin',
        levels=N_iz_lvls, tp=False, tpAlpha=0.2,
        extend='neither', show=False, ext='.pdf',
        params=params)



