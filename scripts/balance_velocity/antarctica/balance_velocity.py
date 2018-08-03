from cslvr import *

mesh_H = 10

# set plot directory :
plt_dir = '../../../images/balance_velocity/antarctica/'

# load the data :
f = HDF5File(mpi_comm_world(), 'dump/vars/state.h5', 'r')

# the balance velocity uses a 2D-model :
model = D2Model(f, out_dir = 'results/')
    
# set the calculated subdomains :
model.set_subdomains(f)

# use the projection of the dataset 'bedmap1' for plotting :
bm1  = DataFactory.get_bedmap1()

model.init_S(f)
model.init_B(f)
model.init_adot(f)
model.init_mask(f)
model.init_U_ob(f,f)
model.init_U_mask(f)

# containers for the direction of flow :
d_U_ob_S_x = Function(model.Q)
d_U_ob_S_y = Function(model.Q)
d_gS_m_U_x = Function(model.Q)
d_gS_m_U_y = Function(model.Q)

# calculate the down-surface gradient :
dSdx = project(-model.S.dx(0))
dSdy = project(-model.S.dx(1))

# convert to numpy arrays :
dSdx_v       = dSdx.vector().array()
dSdy_v       = dSdy.vector().array()
d_U_ob_S_x_v = d_U_ob_S_x.vector().array()
d_U_ob_S_y_v = d_U_ob_S_y.vector().array()
d_gS_m_U_x_v = d_gS_m_U_x.vector().array()
d_gS_m_U_y_v = d_gS_m_U_y.vector().array()
u_ob_v       = model.u_ob.vector().array()
v_ob_v       = model.v_ob.vector().array()

# grounded: down grad(S) -- floating: U_ob  :
d_U_ob_S_x_v[model.shf_dofs] = u_ob_v[model.shf_dofs]
d_U_ob_S_y_v[model.shf_dofs] = v_ob_v[model.shf_dofs]
d_U_ob_S_x_v[model.gnd_dofs] = dSdx_v[model.gnd_dofs]
d_U_ob_S_y_v[model.gnd_dofs] = dSdy_v[model.gnd_dofs]

# everywhere with U observations: U_ob -- everywhere without: down grad(S)  :
d_gS_m_U_x_v[model.Uob_dofs]         = u_ob_v[model.Uob_dofs]
d_gS_m_U_y_v[model.Uob_dofs]         = v_ob_v[model.Uob_dofs]
d_gS_m_U_x_v[model.Uob_missing_dofs] = dSdx_v[model.Uob_missing_dofs]
d_gS_m_U_y_v[model.Uob_missing_dofs] = dSdy_v[model.Uob_missing_dofs]

# assign the numpy arrays back to the containers :
model.assign_variable(d_U_ob_S_x, d_U_ob_S_x_v)
model.assign_variable(d_U_ob_S_y, d_U_ob_S_y_v)
model.assign_variable(d_gS_m_U_x, d_gS_m_U_x_v)
model.assign_variable(d_gS_m_U_y, d_gS_m_U_y_v)

# the imposed direction of flow :
#d = (d_U_ob_S_x, d_U_ob_S_y)
#d = (d_gS_m_U_x, d_gS_m_U_y)
#d = (model.u_ob, model.v_ob)
d = (-model.S.dx(0), -model.S.dx(1))

## plot the observed surface speed :
#U_max  = model.U_ob.vector().max()
#U_min  = model.U_ob.vector().min()
#U_lvls = array([U_min, 2, 10, 20, 50, 100, 200, 500, 1000, U_max])
#plotIce(bm1, model.U_ob, name='U_ob', direc=plt_dir,
#       title=r'$\Vert \mathbf{u}_{ob} \Vert$', cmap='viridis',
#       show=False, levels=U_lvls, tp=False, cb_format='%.1e')

kappas  = [0,5,10]
methods = ['SUPG', 'SSM', 'GLS']

for kappa in kappas:

  for method in methods:

    bv = BalanceVelocity(model, kappa=kappa, stabilization_method=method)
    bv.solve_direction_of_flow(d)
    bv.solve()

    U_max  = model.Ubar.vector().max()
    U_min  = model.Ubar.vector().min()
    U_lvls = array([U_min, 2, 10, 20, 50, 100, 200, 500, 1000, U_max])
    
    name = 'Ubar_%iH_kappa_%i_%s' % (mesh_H, kappa, method)
    tit  = r'$\bar{u}_{%i}$' % kappa
    plotIce(bm1, model.Ubar, name=name, direc=plt_dir,
           title=tit, cmap='viridis',
           show=False, levels=U_lvls, tp=False, cb_format='%.1e')
   
    # calculate the misfit
    misfit = Function(model.Q)
    Ubar_v = model.Ubar.vector().array()
    U_ob_v = model.U_ob.vector().array()
    m_v    = U_ob_v - Ubar_v
    model.assign_variable(misfit, m_v)
   
    m_max  = misfit.vector().max()
    m_min  = misfit.vector().min()
    #m_lvls = array([m_min, -5e2, -1e2, -1e1, -1, 1, 1e1, 1e2, 5e2, m_max])
    m_lvls = array([m_min, -50, -10, -5, -1, 1, 5, 10, 50, m_max])
     
    name = 'misfit_%iH_kappa_%i_%s' % (mesh_H, kappa, method)
    tit  = r'$M_{%i}$' % kappa
    plotIce(bm1, misfit, name=name, direc=plt_dir,
           title=tit, cmap='RdGy',
           show=False, levels=m_lvls, tp=False, cb_format='%.1e')



