from cslvr            import *
from scipy            import random
from fenics           import *
from dolfin_adjoint   import *
import sys


# set the relavent directories :
var_dir = 'dump/vars_jakobshavn_small/'  # directory from gen_vars.py
out_dir = 'dump/jakob_small/inversion_Wc_0.03/'

# create HDF5 files for saving and loading data :
fmeshes = HDF5File(mpi_comm_world(), var_dir + 'submeshes.h5', 'r')
fdata   = HDF5File(mpi_comm_world(), var_dir + 'state.h5',     'r')

# create 3D model for stokes solves :
d3model = D3Model(fdata, out_dir)

# init subdomains and boundary meshes :
d3model.set_subdomains(fdata)
d3model.set_srf_mesh(fmeshes)
d3model.set_bed_mesh(fmeshes)
d3model.set_dvd_mesh(fmeshes)

# initialize the 3D model vars :
d3model.init_S(fdata)                          # surface topography
d3model.init_B(fdata)                          # basal topography
d3model.init_mask(fdata)                       # ice-shelf mask (here no shelf}
d3model.init_q_geo(d3model.ghf)                # geothermal heat flux
d3model.init_T_surface(fdata)                  # surface temperature
d3model.init_adot(fdata)                       # accumulation/ablation
d3model.init_U_ob(fdata, fdata)                # observed velocity
d3model.init_U_mask(fdata)                     # where U_ob data is missing
d3model.init_W(0.0)                            # no water content initially
d3model.init_Wc(0.03)                          # max allowed basal water content
d3model.init_T(d3model.T_surface)              # initial temperature
d3model.init_k_0(1e-3)                         # non-advective enthalpy flux
d3model.solve_hydrostatic_pressure()           # init to hydrostatic pressure
d3model.form_energy_dependent_rate_factor()    # we are thermo-mechanical

# this commented-out section is how we can initialize the data to a previous
# state :

#frstrt = HDF5File(mpi_comm_world(), out_dir + '01/inverted.h5', 'r')
#d3model.init_T(frstrt)
#d3model.init_W(frstrt)
#d3model.init_Fb(frstrt)
#d3model.init_alpha(frstrt)
#d3model.init_U(frstrt)
#d3model.init_p(frstrt)
#d3model.init_theta(frstrt)

# gaps in surface velocity data will create problems for our SIA-approximated
# initial traction values \beta_{SIA}.  To correct this, we use instead of 
# the surface velocity, the balance velocity, in these regions.  See chapter
# on data assimilation for more info on this field 
#===============================================================================

# create a 2D model for balance-velocity :
bedmodel = D2Model(d3model.bedmesh, out_dir)

# initialize the field we need from the 3D mesh :
bedmodel.assign_submesh_variable(bedmodel.S,      d3model.S)
bedmodel.assign_submesh_variable(bedmodel.B,      d3model.B)
bedmodel.assign_submesh_variable(bedmodel.adot,   d3model.adot)

# solve the balance velocity (see the appropriate Chapter) :
bv = BalanceVelocity(bedmodel, kappa=5.0, stabilization_method = 'GLS')
bv . solve_direction_of_flow( [bedmodel.S.dx(0), bedmodel.S.dx(1)] )
bv.solve(annotate=False)

# assign the balance velocity to the 3D model's bed :
d3model.assign_submesh_variable(d3model.d_x,  bedmodel.d_x)
d3model.assign_submesh_variable(d3model.d_y,  bedmodel.d_y)
d3model.assign_submesh_variable(d3model.Ubar, bedmodel.Ubar)

# extrude the bed values up the ice column : 
d_x_e  = d3model.vert_extrude(d3model.d_x,  d='up')
d_y_e  = d3model.vert_extrude(d3model.d_y,  d='up')
Ubar_e = d3model.vert_extrude(d3model.Ubar, d='up')

# set the appropriate variable to be the function extruded :
d3model.init_d_x(d_x_e)
d3model.init_d_y(d_y_e)
d3model.init_Ubar(Ubar_e)

# generate initial traction field from the SIA approximation :
d3model.init_beta_SIA()

# the assimilation process uses the first-order model :
mom    = MomentumDukowiczBP(d3model, linear=False)

# the TMC process uses full-Stokes :
momTMC = MomentumDukowiczStokes(d3model, linear=False)

# the energy-balance physics :
nrg    = Enthalpy(d3model, momTMC, transient=False, use_lat_bc=True)

# this commented out section allows you to restart a simulation from a previous
# momentum optimization :

#frstrt = HDF5File(mpi_comm_world(), out_dir + '02/u_opt.h5', 'r')
#d3model.set_out_dir(out_dir + '02/')
#d3model.init_U(frstrt)
#d3model.init_beta(frstrt)

# thermo-solve callback function is called after every TMC iteration :
def tmc_cb_ftn():
  nrg.calc_PE()#avg=True)
  nrg.calc_vert_avg_strain_heat()
  nrg.calc_vert_avg_W()
  nrg.calc_temp_rat()

# post-adjoint-iteration callback function is called after the TMC data 
# assimilation process is finished :
def adj_post_cb_ftn():
  # solve for optimal vertical velocity :
  mom.solve_vert_velocity(annotate=False)

# after every completed adjoining, save the state of these functions :
adj_save_vars = [d3model.T,
                 d3model.W,
                 d3model.Fb,
                 d3model.Mb,
                 d3model.alpha,
                 d3model.alpha_int,
                 d3model.PE,
                 d3model.Wbar,
                 d3model.Qbar,
                 d3model.temp_rat,
                 d3model.U3,
                 d3model.p,
                 d3model.beta,
                 d3model.theta]

# save these variables after the momentum optimization :
u_opt_save_vars = [d3model.beta, d3model.U3]

# save these variables after the energy optimization :
w_opt_save_vars = [d3model.Fb,   d3model.theta]

# form the momentum cost functional :
mom.form_obj_ftn(integral=d3model.GAMMA_U_GND, kind='log_L2_hybrid', 
                 g1=0.01, g2=5000)

# form the traction-regularization functional :
mom.form_reg_ftn(d3model.beta, integral=d3model.GAMMA_B_GND,
                 kind='TV_Tik_hybrid', alpha_tik=1e-1, alpha_tv=10.0)
#mom.form_reg_ftn(d3model.beta, integral=d3model.GAMMA_B_GND,
#                 kind='TV', alpha=10.0)
#mom.form_reg_ftn(d3model.beta, integral=d3model.GAMMA_B_GND,
#                  kind='Tikhonov', alpha=1e-6)

# form the objective functional for water-flux optimization :
nrg.form_cost_ftn(kind='L2')

# keyword arguments to the energy-optimization function :
wop_kwargs = {'max_iter'            : 350, 
              'bounds'              : (0.0, 100.0),
              'method'              : 'ipopt',
              'adj_save_vars'       : w_opt_save_vars,
              'adj_callback'        : None}
                                    
# keyword arguments to the TMC function :
tmc_kwargs = {'momentum'            : momTMC,
              'energy'              : nrg,
              'wop_kwargs'          : wop_kwargs,
              'callback'            : tmc_cb_ftn,
              'atol'                : 1e2,
              'rtol'                : 1e0,
              'max_iter'            : 5,
              'iter_save_vars'      : None,
              'post_tmc_save_vars'  : None,
              'starting_i'          : 1}

# keyword arguments to the momentum optimization function :
uop_kwargs = {'control'             : d3model.beta,
              'bounds'              : (1e-5, 1e7),
              'method'              : 'ipopt',
              'max_iter'            : 1000,
              'adj_save_vars'       : u_opt_save_vars,
              'adj_callback'        : None,
              'post_adj_callback'   : adj_post_cb_ftn}

# keyword arguments for the TMC data assimilation function : 
ass_kwargs = {'momentum'            : mom,
              'beta_i'              : d3model.beta.copy(True),
              'max_iter'            : 10,
              'tmc_kwargs'          : tmc_kwargs,
              'uop_kwargs'          : uop_kwargs,
              'atol'                : 1.0,
              'rtol'                : 1e-4,
              'initialize'          : True,
              'incomplete'          : True,
              'post_iter_save_vars' : adj_save_vars,
              'post_ini_callback'   : None,
              'starting_i'          : 1}

# assimilate ! :
d3model.assimilate_U_ob(**ass_kwargs) 

# or simply optimize the traction coefficient without thermo :
#mom.optimize_U_ob(**uop_kwargs)
                                    
# or only thermo_solve :
#d3model.thermo_solve(**tmc_kwargs)

 

