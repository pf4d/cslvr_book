from cslvr import *

#e_mode  = 'zero_energy'
e_mode  = 'Fb'
out_dir = 'ps_results/' + e_mode + '/'
plt_dir = '../../images/tmc/plane_strain/' + e_mode + '/'

f      = HDF5File(mpi_comm_world(), out_dir + 'tmc.h5',   'r')
fstate = HDF5File(mpi_comm_world(), out_dir + 'state.h5', 'r')

model  = LatModel(fstate, out_dir = out_dir, use_periodic = False)
model.set_subdomains(fstate)

model.init_T(f)
model.init_W(f)
model.init_Fb(f)
model.init_Mb(f)
model.init_alpha(f)
model.init_alpha_int(f)
model.init_PE(f)
model.init_Wbar(f)
model.init_Qbar(f)
model.init_temp_rat(f)
model.init_U(f)
model.init_p(f)
model.init_beta(f)
model.init_theta(f)

#===============================================================================
# plotting :

figsize = (10,2.2)

U_min  = model.U_mag.vector().min()
U_max  = model.U_mag.vector().max()
U_lvls = array([U_min, 100, 200, 400, 600, U_max])
plot_variable(u = model.U_mag, name = 'U_mag', direc = plt_dir,
              figsize             = figsize,
              title               = r'$\Vert \mathbf{u} \Vert$',
              cmap                = 'viridis',
              levels              = U_lvls,
              show                = False,
              ylabel              = r'$z$',
              equal_axes          = False,
              extend              = 'both',
              cb_format           = '%.1e')

p_min  = model.p.vector().min()
p_max  = model.p.vector().max()
p_lvls = array([p_min, 1e6, 5e6, 1e7, 1.5e7, 2e7, 2.5e7, p_max])
plot_variable(u = model.p, name = 'p', direc = plt_dir,
              figsize             = figsize,
              title               = r'$p$',
              cmap                = 'viridis',
              levels              = p_lvls,
              show                = False,
              ylabel              = r'$z$',
              equal_axes          = False,
              extend              = 'both',
              cb_format           = '%.1e')

beta_lvls = array([0, 200, 400, 600, 800, 1000])
plot_variable(u = model.beta, name = 'beta', direc = plt_dir,
              figsize             = (6,2),
              title               = r'$\beta$',
              cmap                = 'viridis',
              levels              = beta_lvls,
              show                = False,
              ylabel              = r'$z$',
              equal_axes          = False,
              cb_format           = '%g')

T_min  = model.T.vector().min()
T_max  = model.T.vector().max()
T_lvls = array([T_min, 230, 240, 250, 260, 265, 270, T_max])
plot_variable(u = model.T, name = 'T', direc = plt_dir,
              figsize             = figsize,
              title               = r'$T$',
              cmap                = 'viridis',
              levels              = T_lvls,
              show                = False,
              ylabel              = r'$z$',
              equal_axes          = False,
              extend              = 'both',
              cb_format           = '%.1f')

W_min  = model.W.vector().min()
W_max  = model.W.vector().max()
W_lvls = array([0.0, 1e-2, 5e-2, W_max])
plot_variable(u = model.W, name = 'W', direc = plt_dir,
              figsize             = figsize,
              title               = r'$W$',
              cmap                = 'viridis',
              levels              = W_lvls,
              show                = False,
              ylabel              = r'$z$',
              equal_axes          = False,
              extend              = 'both',
              cb_format           = '%.1e')

plot_variable(u = model.theta, name = 'theta', direc = plt_dir,
              figsize             = figsize,
              title               = r'$\theta$',
              cmap                = 'viridis',
              levels              = None,
              show                = False,
              ylabel              = r'$z$',
              equal_axes          = False,
              extend              = 'both',
              cb_format           = '%g')

Mb_min  = model.Mb.vector().min()
Mb_max  = model.Mb.vector().max()
Mb_lvls = array([0.0, 0.2, 0.3, 0.4, 0.5, Mb_max])
plot_variable(u = model.Mb, name = 'Mb', direc = plt_dir,
              figsize             = figsize,
              title               = r'$M_b$',
              cmap                = 'viridis',#'gist_yarg',
              levels              = Mb_lvls,
              show                = False,
              ylabel              = r'$z$',
              equal_axes          = False,
              extend              = 'both',
              cb_format           = '%.1f')

Fb_min  = model.Fb.vector().min()
Fb_max  = model.Fb.vector().max()
Fb_lvls = array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, Fb_max])
plot_variable(u = model.Fb, name = 'Fb', direc = plt_dir,
              figsize             = figsize,
              title               = r'$F_b$',
              cmap                = 'viridis',
              levels              = Fb_lvls,
              show                = False,
              ylabel              = r'$z$',
              equal_axes          = False,
              extend              = 'both',
              cb_format           = '%.1f')



