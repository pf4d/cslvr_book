from cslvr  import *
from scipy  import random

#reg_typ = 'TV'
reg_typ = 'Tikhonov'
out_dir = './L_curve_results_' + reg_typ + '/'
plt_dir = '../../images/data_assimilation/ISMIP_HOM_C/'

a     = 0.5 * pi / 180
L     = 20000
bmax  = 1000

p1    = Point(0.0, 0.0, 0.0)
p2    = Point(L,   L,   1)
mesh  = BoxMesh(p1, p2, 15, 15, 5)

model = D3Model(mesh, out_dir = out_dir, use_periodic = True)

surface = Expression('- x[0] * tan(a)', a=a, 
                     element=model.Q.ufl_element())
bed     = Expression('- x[0] * tan(a) - 1000.0', a=a, 
                     element=model.Q.ufl_element())
beta    = Expression('bmax/2 + bmax/2 * sin(2*pi*x[0]/L) * sin(2*pi*x[1]/L)',
                     bmax=bmax, L=L, element=model.Q.ufl_element())

# calculate the boundaries for integration :
model.calculate_boundaries()

# deform the mesh to the desired geometry :
model.deform_mesh_to_geometry(surface, bed)

# initialize important variables :
model.init_beta(beta)                          # traction
model.init_A(1e-16)                            # isothermal rate-factor

mom = MomentumDukowiczBP(model)
mom.solve(annotate=False)

# add noise with a signal-to-noise ratio of 100 :
snr   = 100.0
u     = Function(model.Q)
v     = Function(model.Q)
assign(u, model.U3.sub(0))
assign(v, model.U3.sub(1))
u_o   = u.vector().array()
v_o   = v.vector().array()
n     = len(u_o)
sig   = model.get_norm(as_vector([u, v]), 'linf')[1] / snr
print_min_max(snr, 'SNR')
print_min_max(sig, 'sigma')
  
u_error = sig * random.randn(n)
v_error = sig * random.randn(n)
u_ob    = u_o + u_error
v_ob    = v_o + v_error

# init the 'observed' velocity :
model.init_U_ob(u_ob, v_ob)
u_ob_ex = model.vert_extrude(model.u_ob, 'down')
v_ob_ex = model.vert_extrude(model.v_ob, 'down')
model.init_U_ob(u_ob_ex, v_ob_ex)

# init the traction to the SIA approximation :
model.init_beta_SIA()

# form the cost functional :
mom.form_obj_ftn(integral=model.GAMMA_U_GND, kind='log_L2_hybrid', 
                 g1=1, g2=1e5)

# solving the incomplete adjoint is more efficient :
mom.linearize_viscosity()

# optimize for beta :
adj_kwargs = {'control'           : model.beta,
              'bounds'            : (1e-5, 1e7),
              'method'            : 'ipopt',
              'max_iter'          : 100,
              'adj_save_vars'     : None,
              'adj_callback'      : None,
              'post_adj_callback' : None}

# regularization parameters :
alphas = [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4]

model.L_curve(alphas        = alphas,
              physics       = mom,
              control       = model.beta,
              int_domain    = model.GAMMA_B_GND,
              adj_ftn       = mom.optimize_U_ob,
              adj_kwargs    = adj_kwargs,
              reg_kind      = reg_typ,
              pre_callback  = None,
              post_callback = None,
              itr_save_vars = None)


