from cslvr import *

out_dir = '../../images/stress_balance/BP/'
in_dir  = './results/'

f  = HDF5File(mpi_comm_world(), in_dir + 'BP.h5', 'r')

model = D3Model(f, out_dir = out_dir)

model.set_subdomains(f)

model.init_U(f)
model.init_p(f)
model.init_N_ii(f)
model.init_N_ij(f)
model.init_N_iz(f)
model.init_N_ji(f)
model.init_N_jj(f)
model.init_N_jz(f)
model.init_N_zi(f)
model.init_N_zj(f)
model.init_N_zz(f)
model.init_M_ii(f)
model.init_M_ij(f)
model.init_M_iz(f)
model.init_M_ji(f)
model.init_M_jj(f)
model.init_M_jz(f)
model.init_M_zi(f)
model.init_M_zj(f)
model.init_M_zz(f)

#===============================================================================
# plotting :

# create the bed and surface meshes :
model.form_bed_mesh()
model.form_srf_mesh()

# create 2D models :
bedmodel = D2Model(model.bedmesh, out_dir)
srfmodel = D2Model(model.srfmesh, out_dir)

srfmodel.assign_submesh_variable(srfmodel.U_mag,  model.U_mag)
bedmodel.assign_submesh_variable(bedmodel.p,      model.p)
srfmodel.assign_submesh_variable(srfmodel.N_ii,   model.N_ii)
srfmodel.assign_submesh_variable(srfmodel.N_ij,   model.N_ij)
srfmodel.assign_submesh_variable(srfmodel.N_iz,   model.N_iz)
srfmodel.assign_submesh_variable(srfmodel.N_ji,   model.N_ji)
srfmodel.assign_submesh_variable(srfmodel.N_jj,   model.N_jj)
srfmodel.assign_submesh_variable(srfmodel.N_jz,   model.N_jz)
srfmodel.assign_submesh_variable(srfmodel.N_zi,   model.N_zi)
srfmodel.assign_submesh_variable(srfmodel.N_zj,   model.N_zj)
srfmodel.assign_submesh_variable(srfmodel.N_zz,   model.N_zz)
srfmodel.assign_submesh_variable(srfmodel.M_ii,   model.M_ii)
srfmodel.assign_submesh_variable(srfmodel.M_ij,   model.M_ij)
srfmodel.assign_submesh_variable(srfmodel.M_iz,   model.M_iz)
srfmodel.assign_submesh_variable(srfmodel.M_ji,   model.M_ji)
srfmodel.assign_submesh_variable(srfmodel.M_jj,   model.M_jj)
srfmodel.assign_submesh_variable(srfmodel.M_jz,   model.M_jz)
srfmodel.assign_submesh_variable(srfmodel.M_zi,   model.M_zi)
srfmodel.assign_submesh_variable(srfmodel.M_zj,   model.M_zj)
srfmodel.assign_submesh_variable(srfmodel.M_zz,   model.M_zz)

U_min    = srfmodel.U_mag.vector().min()
U_max    = srfmodel.U_mag.vector().max()
         
p_min    = bedmodel.p.vector().min()
p_max    = bedmodel.p.vector().max()

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

U_lvls = array([U_min, 91, 92, 93, 94, 95, 96, 97, 98, 99, U_max])

p_lvls = array([4e6, 5e6, 6e6, 7e6, 8e6, 9e6, 1e7, 1.1e7, 1.2e7, p_max])

# membrane stresses :
N_ii_lvls = np.array([N_ii_min, -2e7, -1e7, -5e6, -1e6, 
                      1e6, 5e6, 1e7, 2e7, N_ii_max])
N_ij_lvls = np.array([N_ij_min, -3e6, -2e6, -1e6, -5e5,
                      5e5, 1e6, 2e6, 3e6, N_ij_max])
N_iz_lvls = np.array([N_iz_min, 1.3e7, 2e7, 3e7, 4e7, 
                      5e7, 6e7, 7e7, 8e7, N_iz_max])

N_jj_lvls = np.array([N_jj_min, -1.5e7, -1e7, -5e6, -1e6, 
                      1e6, 5e6, 1e7, 1.5e7, N_jj_max])
N_jz_lvls = np.array([N_jz_min, -2e7, -1e7, -5e6, -1e6, 
                      1e6, 5e6, 1e7, 2e7, N_jz_max])

N_zz_lvls = np.array([N_zz_min, -4e6, -2e6, -5e5, -1e5, 
                      1e5, 5e5, 2e6, 4e6, N_zz_max])

## membrane stresses :
#N_ii_lvls = np.array([N_ii_min, -2e7, -1e7, -5e6, -1e6, 
#                      1e6, 5e6, 1e7, 2e7, N_ii_max])
#N_ij_lvls = np.array([N_ij_min, -7.5e6, -5e6, -2.5e6, -1e6, 
#                      1e6, 2.5e6, 5e6, 7.5e6, N_ij_max])
#N_iz_lvls = np.array([N_iz_min, 2.1e7, 2.5e7, 3e7, 3.5e7, 4e7, 
#                      4.5e7, 5e7, 5.5e7, N_iz_max])
#
#N_jj_lvls = np.array([N_jj_min, -1.5e7, -1e7, -5e6, -1e6, 
#                      1e6, 5e6, 1e7, 1.5e7, N_jj_max])
#N_jz_lvls = np.array([N_jz_min, -1.5e6, -1e6, -5e5, -1e5,
#                      1e5, 5e5, 1e6, 1.5e6, N_jz_max])
#
#N_zz_lvls = np.array([N_zz_min, -1.25e7, -1e7, -5e6, -1e6, 
#                      1e6, 5e6, 1e7, 1.25e7, N_zz_max])

# membrane stress-balance :
M_ii_lvls = np.array([M_ii_min, -2.5e4, -2e4, -1.5e4, -5e3,
                      5e3, 1.5e4, 2e4, 2.5e4, M_ii_max])
M_ij_lvls = np.array([M_ij_min, -8e3, -5e3, -2.5e3, -5e2,
                      5e2, 2.5e3, 5e3, 8e3, M_ij_max])
M_iz_lvls = np.array([M_iz_min, -6.4e4, -6.2e4, -6.0e4, -5.8e4,
                      -5.6e4, -5.4e4, M_iz_max])

M_ji_lvls = np.array([M_ji_min, -8e3, -5e3, -2.5e3, -5e2,
                      5e2, 2.5e3, 5e3, 8e3, M_ji_max])
M_jj_lvls = np.array([M_jj_min, -1e4, -7.5e3, -5e3, -1e3,
                      1e3, 5e3, 7.5e3, 1e4, M_jj_max])
M_jz_lvls = np.array([M_jz_min, -5e3, -2.5e3, -1e3, -1e2,
                      1e2, 1e3, 2.5e3, 5e3, M_jz_max])

M_zi_lvls = np.array([M_zi_min, -3e4, -2e4, -1e4, -1e3,
                      1e3, 1e4, 2e4, 3e4, M_zi_max])
M_zj_lvls = np.array([M_zj_min, -1.5e3, -1e3, -5e2, -1e2, 
                      1e2, 5e2, 1e3, 1.5e3, M_zj_max])
M_zz_lvls = np.array([M_zz_min, -2e4, -1.5e4, -1e4, -5e3, 
                      5e3, 1e4, 1.5e4, 2e4, M_zz_max])


#===============================================================================
# plot :

plot_variable(u                   = srfmodel.U_mag,
              name                = 'U_mag',
              levels              = U_lvls,
              direc               = out_dir,
              cmap                = 'viridis',
              tp                  = True,
              show                = False,
              extend              = 'neither',
              cb_format           = '%g')

plot_variable(u                   = bedmodel.p,
              name                = 'p',
              levels              = p_lvls,
              direc               = out_dir,
              cmap                = 'viridis',
              tp                  = True,
              show                = False,
              extend              = 'min',
              cb_format           = '%.1e')

plot_variable(u                   = srfmodel.N_ii,
              name                = 'N_ii',
              levels              = N_ii_lvls,
              direc               = out_dir,
              cmap                = 'RdGy',
              tp                  = True,
              show                = False,
              extend              = 'neither',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

plot_variable(u                   = srfmodel.N_ij,
              name                = 'N_ij',
              levels              = N_ij_lvls,
              direc               = out_dir,
              cmap                = 'RdGy',
              tp                  = True,
              show                = False,
              extend              = 'neither',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

plot_variable(u                   = srfmodel.N_iz,
              name                = 'N_iz',
              levels              = N_iz_lvls,
              direc               = out_dir,
              cmap                = 'gist_yarg',
              tp                  = True,
              show                = False,
              extend              = 'neither',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

plot_variable(u                   = srfmodel.N_ji,
              name                = 'N_ji',
              levels              = N_ij_lvls,
              direc               = out_dir,
              cmap                = 'RdGy',
              tp                  = True,
              show                = False,
              extend              = 'neither',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

plot_variable(u                   = srfmodel.N_jj,
              name                = 'N_jj',
              levels              = N_jj_lvls,
              direc               = out_dir,
              cmap                = 'RdGy',
              tp                  = True,
              show                = False,
              extend              = 'neither',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

plot_variable(u                   = srfmodel.N_jz,
              name                = 'N_jz',
              levels              = N_jz_lvls,
              direc               = out_dir,
              cmap                = 'RdGy',
              tp                  = True,
              show                = False,
              extend              = 'neither',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = srfmodel.N_zi,
              name                = 'N_zi',
              levels              = N_iz_lvls,
              direc               = out_dir,
              cmap                = 'gist_yarg',
              tp                  = True,
              show                = False,
              extend              = 'neither',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = srfmodel.N_zj,
              name                = 'N_zj',
              levels              = N_jz_lvls,
              direc               = out_dir,
              cmap                = 'RdGy',
              tp                  = True,
              show                = False,
              extend              = 'neither',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = srfmodel.N_zz,
              name                = 'N_zz',
              levels              = N_zz_lvls,
              direc               = out_dir,
              cmap                = 'RdGy',
              tp                  = True,
              show                = False,
              extend              = 'neither',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
                                  
                                  
plot_variable(u                   = srfmodel.M_ii,
              name                = 'M_ii',
              levels              = M_ii_lvls,
              direc               = out_dir,
              cmap                = 'RdGy',
              tp                  = True,
              show                = False,
              extend              = 'neither',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = srfmodel.M_ij,
              name                = 'M_ij',
              levels              = M_ij_lvls,
              direc               = out_dir,
              cmap                = 'RdGy',
              tp                  = True,
              show                = False,
              extend              = 'neither',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = srfmodel.M_iz,
              name                = 'M_iz',
              levels              = M_iz_lvls,
              direc               = out_dir,
              cmap                = 'Reds',
              tp                  = True,
              show                = False,
              extend              = 'neither',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = srfmodel.M_ji,
              name                = 'M_ji',
              levels              = M_ji_lvls,
              direc               = out_dir,
              cmap                = 'RdGy',
              tp                  = True,
              show                = False,
              extend              = 'neither',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = srfmodel.M_jj,
              name                = 'M_jj',
              levels              = M_jj_lvls,
              direc               = out_dir,
              cmap                = 'RdGy',
              tp                  = True,
              show                = False,
              extend              = 'neither',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = srfmodel.M_jz,
              name                = 'M_jz',
              levels              = M_jz_lvls,
              direc               = out_dir,
              cmap                = 'RdGy',
              tp                  = True,
              show                = False,
              extend              = 'neither',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = srfmodel.M_zi,
              name                = 'M_zi',
              levels              = M_zi_lvls,
              direc               = out_dir,
              cmap                = 'RdGy',
              tp                  = True,
              show                = False,
              extend              = 'neither',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = srfmodel.M_zj,
              name                = 'M_zj',
              levels              = M_zj_lvls,
              direc               = out_dir,
              cmap                = 'RdGy',
              tp                  = True,
              show                = False,
              extend              = 'neither',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = srfmodel.M_zz,
              name                = 'M_zz',
              levels              = M_zz_lvls,
              direc               = out_dir,
              cmap                = 'RdGy',
              tp                  = True,
              show                = False,
              extend              = 'neither',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')



