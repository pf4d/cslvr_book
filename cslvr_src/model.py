from fenics               import *
from dolfin_adjoint       import *
from cslvr.io             import print_text, get_text, print_min_max
from copy                 import copy
import numpy              as np
import matplotlib.pyplot  as plt
import matplotlib         as mpl
import sys
import os
import re


class Model(object):
  """ 
  The basic model from which each of these inherit :
   * :class:`~latmodel.LatModel`       - plane strain model
   * :class:`~hybridmodel.HybridModel` - L1L2 "hybrid" model
   * :class:`~d1model.D1Model`         - 1D firn model
   * :class:`~d2model.D2Model`         - 2D model (SSA, SIA, balance velocity)
   * :class:`~d3model.D3Model`         - 3D model (first-order, full-Stokes)

  Args:
    :mesh:         the :class:`~fenics.Mesh` instance.
                         
    :out_dir:      string location for the output directory, default
                   ``./results/``.

    :use_periodic: boolean, use periodic boundaries or not, default 
                   ``False``.
  """

  OMEGA_GND   = 0   # internal cells over bedrock
  OMEGA_FLT   = 1   # internal cells over water
  GAMMA_S_GND = 2   # grounded upper surface
  GAMMA_B_GND = 3   # grounded lower surface (bedrock)
  GAMMA_S_FLT = 6   # shelf upper surface
  GAMMA_B_FLT = 5   # shelf lower surface
  GAMMA_L_DVD = 7   # basin divides
  GAMMA_L_OVR = 4   # terminus over water
  GAMMA_L_UDR = 10  # terminus under water
  GAMMA_U_GND = 8   # grounded surface with U observations
  GAMMA_U_FLT = 9   # shelf surface with U observations
  
  # external boundaries :
  ext_boundaries = {GAMMA_S_GND : 'grounded upper surface',
                    GAMMA_B_GND : 'grounded lower surface (bedrock)',
                    GAMMA_S_FLT : 'shelf upper surface',
                    GAMMA_B_FLT : 'shelf lower surface',
                    GAMMA_L_DVD : 'basin divides',
                    GAMMA_L_OVR : 'terminus over water',
                    GAMMA_L_UDR : 'terminus under water',
                    GAMMA_U_GND : 'grounded upper surface with U observations',
                    GAMMA_U_FLT : 'shelf upper surface with U observations'}

  # internal boundaries :
  int_boundaries = {OMEGA_GND   : 'internal cells laying over bedrock',
                    OMEGA_FLT   : 'internal cells laying over water'}

  # union :
  boundaries = dict(ext_boundaries, **int_boundaries)
  
  def __init__(self, mesh, out_dir='./results/', order=1, use_periodic=False):
    """
    Create and instance of the model.
    """
    self.this = super(type(self), self)  # pointer to this base class
  
    s = "::: INITIALIZING BASE MODEL :::"
    print_text(s, cls=self.this)
    
    parameters['form_compiler']['quadrature_degree']  = 2
    parameters["std_out_all_processes"]               = False
    parameters['form_compiler']['cpp_optimize']       = True

    PETScOptions.set("mat_mumps_icntl_14", 100.0)

    self.order       = order
    self.out_dir     = out_dir
    self.MPI_rank    = MPI.rank(mpi_comm_world())
    self.use_periodic_boundaries = use_periodic
    
    self.generate_constants()
    self.set_mesh(mesh)
    self.generate_function_spaces(order, use_periodic)
    self.initialize_variables()

  def color(self):
    return '148'

  def generate_constants(self):
    """
    Initializes important constants.
    """
    s = "::: generating constants :::"
    print_text(s, cls=self.this)

    spy = 31556926.0
    ghf = 0.042 * spy  # W/m^2 = J/(s*m^2) = spy * J/(a*m^2)

    # Constants :
    self.kcHh    = Constant(3.7e-9)
    self.kcHh.rename('kcHh', 'creep coefficient high')

    self.kcLw    = Constant(9.2e-9)
    self.kcLw.rename('kcLw', 'creep coefficient low ')

    self.kg      = Constant(1.3e-7)
    self.kg.rename('kg', 'grain growth coefficient')

    self.Ec      = Constant(6e4)
    self.Ec.rename('Ec', 'act. energy for water in ice')

    self.Eg      = Constant(42.4e3)
    self.Eg.rename('Eg', 'act. energy for grain growth')

    self.thetasp = Constant(2009.0 * 273.15)
    self.thetasp.rename('thetasp', 'Internal energy of ice at Tw')
    
    self.etaw    = Constant(1.787e-3)
    self.etaw.rename('etaw', 'Dynamic viscosity of water at Tw')

    self.r       = Constant(0.0)
    self.r.rename('r', 'thickness exponent in sliding law')

    self.eps_reg = Constant(1e-15)
    self.eps_reg.rename('eps_reg', 'strain rate regularization parameter')

    self.n       = Constant(3.0)
    self.n.rename('n', 'viscosity nonlinearity parameter')

    self.spy     = Constant(spy)
    self.spy.rename('spy', 'seconds per year')

    self.rhoi    = Constant(910.0)
    self.rhoi.rename('rhoi', 'ice density')

    self.rhow    = Constant(1000.0)
    self.rhow.rename('rhow', 'water density')

    self.rhosw   = Constant(1028.0)
    self.rhosw.rename('rhosw', 'sea-water density')
    
    self.rhom    = Constant(550.0)
    self.rhom.rename('rhom', 'firn pore close-off density')

    self.rhoc    = Constant(815.0)
    self.rhoc.rename('rhoc', 'firn density critical value')

    self.g       = Constant(9.80665)
    self.g.rename('g', 'gravitational acceleration')

    self.a0      = Constant(5.45e10)
    self.a0.rename('a0', 'ice hardness limit')

    self.Q0      = Constant(13.9e4)
    self.Q0.rename('Q0', 'ice activation energy')

    self.R       = Constant(8.3144621)
    self.R.rename('R', 'universal gas constant')

    self.ki      = Constant(2.1)
    self.ki.rename('ki', 'thermal conductivity of ice')

    self.kw      = Constant(0.561)
    self.kw.rename('kw', 'thermal conductivity of water')

    self.ci      = Constant(2009.0)
    self.ci.rename('ci', 'heat capacity of ice')
    
    self.cw      = Constant(4217.6)
    self.cw.rename('cw', 'Heat capacity of water at Tw')

    self.L       = Constant(3.3355e5)
    self.L.rename('L', 'latent heat of ice')

    self.ghf     = Constant(ghf)
    self.ghf.rename('ghf', 'geothermal heat flux')

    self.gamma   = Constant(9.8e-8)
    self.gamma.rename('gamma', 'pressure melting point depth dependence')

    self.nu      = Constant(3.5e3)
    self.nu.rename('nu', 'moisture diffusivity')

    self.T_w     = Constant(273.15)
    self.T_w.rename('T_w', 'Triple point of water')

    self.a_T_l   = Constant(3.985e-13*spy)
    self.a_T_l.rename('a_T_l', 'lower bound of flow-rate constant')

    self.a_T_u   = Constant(1.916e3*spy)
    self.a_T_u.rename('a_T_u', 'upper bound of flow-rate constant')

    self.Q_T_l   = Constant(6e4)
    self.Q_T_l.rename('Q_T_l', 'lower bound of ice activation energy')

    self.Q_T_u   = Constant(13.9e4)
    self.Q_T_u.rename('Q_T_u', 'upper bound of ice activation energy')
  
  def set_measures(self, ff, cf):
    """
    set the new measure space for facets (self.ds) and cells (self.dx) for
    the boundaries marked by FacetFunction *ff* and CellFunction *cf*.
    """
    # calculate the number of cells and facets that are of a certain type
    # for determining Dirichlet boundaries :
    self.N_OMEGA_GND   = sum(self.cf.array() == self.OMEGA_GND)
    self.N_OMEGA_FLT   = sum(self.cf.array() == self.OMEGA_FLT)
    self.N_GAMMA_S_GND = sum(self.ff.array() == self.GAMMA_S_GND)
    self.N_GAMMA_B_GND = sum(self.ff.array() == self.GAMMA_B_GND)
    self.N_GAMMA_S_FLT = sum(self.ff.array() == self.GAMMA_S_FLT)
    self.N_GAMMA_B_FLT = sum(self.ff.array() == self.GAMMA_B_FLT)
    self.N_GAMMA_L_DVD = sum(self.ff.array() == self.GAMMA_L_DVD)
    self.N_GAMMA_L_OVR = sum(self.ff.array() == self.GAMMA_L_OVR)
    self.N_GAMMA_L_UDR = sum(self.ff.array() == self.GAMMA_L_UDR)
    self.N_GAMMA_U_GND = sum(self.ff.array() == self.GAMMA_U_GND)
    self.N_GAMMA_U_FLT = sum(self.ff.array() == self.GAMMA_U_FLT)

    # create new measures of integration :
    self.ds      = Measure('ds', subdomain_data=self.ff)
    self.dx      = Measure('dx', subdomain_data=self.cf)
    
    self.dx_g    = self.dx(0)                # internal above grounded
    self.dx_f    = self.dx(1)                # internal above floating
    self.dBed_g  = self.ds(3)                # grounded bed
    self.dBed_f  = self.ds(5)                # floating bed
    self.dBed    = self.ds(3) + self.ds(5)   # bed
    self.dSrf_gu = self.ds(8)                # grounded with U observations
    self.dSrf_fu = self.ds(9)                # floating with U observations
    self.dSrf_u  = self.ds(8) + self.ds(9)   # surface with U observations
    self.dSrf_g  = self.ds(2) + self.ds(8)   # surface of grounded ice
    self.dSrf_f  = self.ds(6) + self.ds(9)   # surface of floating ice
    self.dSrf    =   self.ds(6) + self.ds(2) \
                   + self.ds(8) + self.ds(9) # surface
    self.dLat_d  = self.ds(7)                # lateral divide
    self.dLat_to = self.ds(4)                # lateral terminus overwater
    self.dLat_tu = self.ds(10)               # lateral terminus underwater
    self.dLat_t  = self.ds(4) + self.ds(10)  # lateral terminus
    self.dLat    =   self.ds(4) + self.ds(7) \
                   + self.ds(10)             # lateral

  def generate_pbc(self):
    """
    return a :class:`fenics.SubDomain` of periodic lateral boundaries.
    """
    raiseNotDefined()
    
  def set_mesh(self, f):
    """
    Sets the ``mesh`` instance to *f*, either a :class:`fenics.Mesh` or  ``.h5``
    file with a mesh saved with name ``mesh``.
    """
    s = "::: setting mesh :::"
    print_text(s, cls=self.this)

    if isinstance(f, dolfin.cpp.io.HDF5File):
      self.mesh = Mesh()
      f.read(self.mesh, 'mesh', False)

    elif isinstance(f, dolfin.cpp.mesh.Mesh):
      self.mesh = f

    self.dim   = self.mesh.ufl_cell().topological_dimension()

  def calculate_boundaries(self):
    """
    Determines the boundaries of the current ``self.mesh``.
    """
    raiseNotDefined()
  
  def set_out_dir(self, out_dir):
    """
    Set the output directory to string *out_dir*.
    """
    self.out_dir = out_dir
    s = "::: output directory changed to '%s' :::" % out_dir
    print_text(s, cls=self.this)

  def generate_function_spaces(self, order=1, use_periodic=False):
    """
    Generates the finite-element function spaces used by all children of this
    :class:`Model`.

    Arg:

      :use_periodic: boolean to use periodic boundaries along lateral boundary.

    """
    s = "::: generating fundamental function spaces of order %i :::" % order
    print_text(s, cls=self.this)

    if use_periodic:
      self.generate_pbc()
    else:
      self.pBC = None
    self.Q      = FunctionSpace(self.mesh,      "CG", order, 
                                constrained_domain=self.pBC)
    Qe          = FiniteElement("CG", self.mesh.ufl_cell(), order)
    self.Q2     = FunctionSpace(self.mesh, MixedElement([Qe]*2),
                                constrained_domain=self.pBC)
    self.Q3     = FunctionSpace(self.mesh, MixedElement([Qe]*3),
                                constrained_domain=self.pBC)
    self.Q4     = FunctionSpace(self.mesh, MixedElement([Qe]*4),
                                constrained_domain=self.pBC)
    self.Q5     = FunctionSpace(self.mesh, MixedElement([Qe]*5),
                                constrained_domain=self.pBC)
    #self.Q2     = MixedFunctionSpace([self.Q]*2)
    #self.Q3     = MixedFunctionSpace([self.Q]*3)
    #self.Q4     = MixedFunctionSpace([self.Q]*4)
    #self.Q5     = MixedFunctionSpace([self.Q]*5)
    self.Q_non_periodic   = FunctionSpace(self.mesh, "CG", order)
    self.Q3_non_periodic  = FunctionSpace(self.mesh, MixedElement([Qe]*3))
    self.V                = VectorFunctionSpace(self.mesh, "CG", order)

    s = "    - fundamental function spaces created - "
    print_text(s, cls=self.this)
  
  def init_S(self, S, cls=None):
    r"""
    Set surface topography :math:`S` to *S*. 
    """
    if cls is None:
      cls = self.this
    s = "::: initializng surface topography :::"
    print_text(s, cls=cls)
    self.assign_variable(self.S, S, cls=cls)

  def init_B(self, B, cls=None):
    r"""
    Set bed topography :math:`B` to *B*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializng bed topography :::"
    print_text(s, cls=cls)
    self.assign_variable(self.B, B, cls=cls)
  
  def init_p(self, p, cls=None):
    r"""
    Set pressure :math:`p` to *p*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing pressure :::"
    print_text(s, cls=cls)
    self.assign_variable(self.p, p, cls=cls)
  
  def init_theta(self, theta, cls=None):
    r"""
    Set internal energy :math:`\theta` to *theta*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing internal energy :::"
    print_text(s, cls=cls)
    self.assign_variable(self.theta, theta, cls=cls)
    # init pressure-melting temperature thing :
    theta_v = self.theta.vector().array()
    T_v     = (-146.3 + np.sqrt(146.3**2 + 2*7.253*theta_v)) / 7.253
    T_v[T_v > self.T_w(0)] = self.T_w(0)
    self.init_Tp(T_v, cls=cls)
  
  def init_theta_app(self, theta_app, cls=None):
    r"""
    Set the internal energy approximation :math:`\theta_{app}` to *theta_app*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing internal energy approximation :::"
    print_text(s, cls=cls)
    self.assign_variable(self.theta_app, theta_app, cls=cls)
  
  def init_theta_surface(self, theta_surface, cls=None):
    r"""
    Set the surface internal energy :math:`\theta_S` to *theta_surface*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing surface energy :::"
    print_text(s, cls=cls)
    self.assign_variable(self.theta_surface, theta_surface, cls=cls)
  
  def init_theta_float(self, theta_float, cls=None):
    r"""
    Set the floating ice internal energy :math:`\theta_{sea}` to *theta_float*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing floating bed energy :::"
    print_text(s, cls=cls)
    self.assign_variable(self.theta_float, theta_float, cls=cls)
  
  def init_T(self, T, cls=None):
    r"""
    Set temperature :math:`T` to *T*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing absolute temperature :::"
    print_text(s, cls=cls)
    self.assign_variable(self.T, T, cls=cls)
  
  def init_Tp(self, Tp, cls=None):
    r"""
    Set pressure-adjusted temperature :math:`T_p = T + \gamma p` to *Tp*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing pressure-adjusted temperature :::"
    print_text(s, cls=cls)
    self.assign_variable(self.Tp, Tp, cls=cls)
  
  def init_W(self, W, cls=None):
    r"""
    Set water content :math:`W` to *W*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing water content :::"
    print_text(s, cls=cls)
    self.assign_variable(self.W, W, cls=cls)
  
  def init_Wc(self, Wc, cls=None):
    r"""
    Set maximum observed water content :math:`W_c` to *Wc*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing maximum water content :::"
    print_text(s, cls=cls)
    self.assign_variable(self.Wc, Wc, cls=cls)
  
  def init_Mb(self, Mb, cls=None):
    r"""
    Set basal melting rate :math:`M_b` to *Mb*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing basal melt rate :::"
    print_text(s, cls=cls)
    self.assign_variable(self.Mb, Mb, cls=cls)
  
  def init_adot(self, adot, cls=None):
    r"""
    Set accumulation/ablation :math:`\dot{a}` to *adot*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing accumulation :::"
    print_text(s, cls=cls)
    self.assign_variable(self.adot, adot, cls=cls)
  
  def init_beta(self, beta, cls=None):
    r"""
    Set basal traction :math:`\beta` to *beta*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing basal traction coefficient :::"
    print_text(s, cls=cls)
    self.assign_variable(self.beta, beta, cls=cls)
  
  def init_A(self, A, cls=None):
    r"""
    Set rate factor :math:`A` to *A*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing rate factor over grounded and shelves :::"
    print_text(s, cls=cls)
    self.assign_variable(self.A, A, cls=cls)
    self.init_A_shf(A, cls=cls)
    self.init_A_gnd(A, cls=cls)
  
  def init_A_shf(self, A_shf, cls=None):
    r"""
    Set overlying shelf rate factor :math:`A |_{sea}` to *A_shf*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing rate factor over shelves :::"
    print_text(s, cls=cls)
    self.assign_variable(self.A_shf, A_shf, cls=cls)
  
  def init_A_gnd(self, A_gnd, cls=None):
    r"""
    Set overlying grounded rate factor :math:`A |_{gnd}` to *A_gnd*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing rate factor over grounded ice :::"
    print_text(s, cls=cls)
    self.assign_variable(self.A_gnd, A_gnd, cls=cls)
    
  def init_E(self, E, cls=None):
    r"""
    Set flow-enhancement factor :math:`E` to *E*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing enhancement factor over grounded and shelves :::"
    print_text(s, cls=cls)
    self.assign_variable(self.E, E, cls=cls)
    self.init_E_shf(E, cls=cls)
    self.init_E_gnd(E, cls=cls)
  
  def init_E_shf(self, E_shf, cls=None):
    r"""
    Set overlying shelf flow-enhancement factor :math:`E |_{shf}` to *E_shf*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing enhancement factor over shelves :::"
    print_text(s, cls=cls)
    self.assign_variable(self.E_shf, E_shf, cls=cls)
  
  def init_E_gnd(self, E_gnd, cls=None):
    r"""
    Set overlying grounded flow-enhancement factor :math:`E |_{gnd}` to *E_gnd*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing enhancement factor over grounded ice :::"
    print_text(s, cls=cls)
    self.assign_variable(self.E_gnd, E_gnd, cls=cls)
  
  def init_eta(self, eta, cls=None):
    r"""
    Set viscosity :math:`\eta` to *eta*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing viscosity :::"
    print_text(s, cls=cls)
    self.assign_variable(self.eta, eta, cls=cls)
  
  def init_etabar(self, etabar, cls=None):
    r"""
    Set vertically averaged viscosity :math:`\bar{\eta}` to *etabar*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing vertically averaged viscosity :::"
    print_text(s, cls=cls)
    self.assign_variable(self.etabar, etabar, cls=cls)
  
  def init_ubar(self, ubar, cls=None):
    r"""
    Set vertically averaged x-component of velocity :math:`\bar{u}` to *ubar*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing vertically averaged x-component of velocity :::"
    print_text(s, cls=cls)
    self.assign_variable(self.ubar, ubar, cls=cls)
  
  def init_vbar(self, vbar, cls=None):
    r"""
    Set vertically averaged y-component of velocity :math:`\bar{v}` to *vbar*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing vertically averaged y-component of velocity :::"
    print_text(s, cls=cls)
    self.assign_variable(self.vbar, vbar, cls=cls)
    
  def init_wbar(self, wbar, cls=None):
    r"""
    Set vertically averaged z-component of velocity :math:`\bar{w}` to *wbar*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing vertically averaged z-component of velocity :::"
    print_text(s, cls=cls)
    self.assign_variable(self.wbar, wbar, cls=cls)
  
  def init_T_surface(self, T_s, cls=None):
    r"""
    Set surface temperature :math:`T_S` to *T_s*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing surface temperature :::"
    print_text(s, cls=cls)
    self.assign_variable(self.T_surface, T_s, cls=cls)
  
  def init_q_geo(self, q_geo, cls=None):
    r"""
    Set geothermal heat flux :math:`q_{geo}` to *q_geo*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing geothermal heat flux :::"
    print_text(s, cls=cls)
    self.assign_variable(self.q_geo, q_geo, cls=cls)
  
  def init_q_fric(self, q_fric, cls=None):
    r"""
    Set traction heat flux :math:`q_{fric}` to *q_fric*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing basal friction heat flux :::"
    print_text(s, cls=cls)
    self.assign_variable(self.q_fric, q_fric, cls=cls)
  
  def init_gradT_B(self, gradT_B, cls=None):
    r"""
    Set basal temperature gradient 
    :math:`\left( k \nabla T \right) \cdot \mathbf{n}` to *gradT_B*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing basal temperature flux :::"
    print_text(s, cls=cls)
    self.assign_variable(self.gradT_B, gradT_B, cls=cls)
  
  def init_gradTm_B(self, gradTm_B, cls=None):
    r"""
    Set basal temperature melting gradient 
    :math:`\left( k \nabla T_m \right) \cdot \mathbf{n}` to *gradTm_B*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing basal temperature-melting flux :::"
    print_text(s, cls=cls)
    self.assign_variable(self.gradTm_B, gradTm_B, cls=cls)
  
  def init_u(self, u, cls=None):
    r"""
    Set x-component of velocity :math:`u` to *u*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing x-component of velocity :::"
    print_text(s, cls=cls)
    u_t = Function(self.Q_non_periodic, name='u_t')
    self.assign_variable(u_t, u, cls=cls)
    self.assx.assign(self.u, u_t, annotate=False)
  
  def init_v(self, v, cls=None):
    r"""
    Set y-component of velocity :math:`v` to *v*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing y-component of velocity :::"
    print_text(s, cls=cls)
    v_t = Function(self.Q_non_periodic, name='v_t')
    self.assign_variable(v_t, v, cls=cls)
    self.assy.assign(self.v, v_t, annotate=False)
  
  def init_w(self, w, cls=None):
    r"""
    Set z-component of velocity :math:`w` to *w*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing z-component of velocity :::"
    print_text(s, cls=cls)
    w_t = Function(self.Q_non_periodic, name='w_t')
    self.assign_variable(w_t, w, cls=cls)
    self.assz.assign(self.w, w_t, annotate=False)

  def init_U(self, U, cls=None):
    r"""
    Set velocity vector :math:`\mathbf{u}` to *U*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing velocity :::"
    print_text(s, cls=cls)
    self.assign_variable(self.U3, U, cls=cls)
    self.init_U_mag(self.U3, cls)

  def init_U_mag(self, U, cls=None):
    r"""
    Set velocity vector magnitude :math:`\Vert \mathbf{u} \Vert` from vector
    *U*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing velocity magnitude :::"
    print_text(s, cls=cls)
    # fenics issue #405 bug workaround :
    if self.use_periodic_boundaries:
      u      = Function(self.Q)
      v      = Function(self.Q)
      w      = Function(self.Q)
      assign(u, U.sub(0))
      assign(v, U.sub(1))
      assign(w, U.sub(2))
    else:
      u,v,w  = U.split(True)
    u_v      = u.vector().array()
    v_v      = v.vector().array()
    w_v      = w.vector().array()
    U_mag_v  = np.sqrt(u_v**2 + v_v**2 + w_v**2 + DOLFIN_EPS)
    self.assign_variable(self.U_mag, U_mag_v, cls=cls)
  
  def init_U_ob(self, u_ob, v_ob, cls=None):
    r"""
    Set horizontal velocity observation vector :math:`\mathbf{u}_{ob}`
    to *U_ob*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing surface velocity :::"
    print_text(s, cls=cls)
    self.assign_variable(self.u_ob, u_ob, cls=cls)
    self.assign_variable(self.v_ob, v_ob, cls=cls)
    u_v      = self.u_ob.vector().array()
    v_v      = self.v_ob.vector().array()
    U_mag_v  = np.sqrt(u_v**2 + v_v**2 + 1e-16)
    self.assign_variable(self.U_ob, U_mag_v, cls=cls)
  
  def init_Ubar(self, Ubar, cls=None):
    r"""
    Set balance velocity :math:`\Vert \bar{\mathbf{u}} \Vert` to *Ubar*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing balance velocity :::"
    print_text(s, cls=cls)
    self.assign_variable(self.Ubar, Ubar, cls=cls)
  
  def init_u_lat(self, u_lat, cls=None):
    r"""
    Set x-component of lateral velocity :math:`u_D` to *u_lat*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing u lateral boundary condition :::"
    print_text(s, cls=cls)
    self.assign_variable(self.u_lat, u_lat, cls=cls)
  
  def init_v_lat(self, v_lat, cls=None):
    r"""
    Set y-component of lateral velocity :math:`v_D` to *v_lat*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing v lateral boundary condition :::"
    print_text(s, cls=cls)
    self.assign_variable(self.v_lat, v_lat, cls=cls)
  
  def init_w_lat(self, w_lat, cls=None):
    r"""
    Set z-component of lateral velocity :math:`w_D` to *w_lat*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing w lateral boundary condition :::"
    print_text(s, cls=cls)
    self.assign_variable(self.w_lat, w_lat, cls=cls)
  
  def init_mask(self, mask, cls=None):
    r"""
    Set shelf mask :math:`M` to *M*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing shelf mask :::"
    print_text(s, cls=cls)
    self.assign_variable(self.mask, mask, cls=cls)
    self.shf_dofs = np.where(self.mask.vector().array() == 2.0)[0]
    self.gnd_dofs = np.where(self.mask.vector().array() == 1.0)[0]
  
  def init_U_mask(self, U_mask, cls=None):
    r"""
    Set velocity observation mask :math:`M_{\mathbf{u}}` to *U_mask*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing velocity mask :::"
    print_text(s, cls=cls)
    self.assign_variable(self.U_mask, U_mask, cls=cls)
    self.Uob_dofs         = np.where(self.U_mask.vector().array() == 1.0)[0]
    self.Uob_missing_dofs = np.where(self.U_mask.vector().array() == 0.0)[0]
  
  def init_lat_mask(self, lat_mask, cls=None):
    r"""
    Set lateral boundary mask :math:`M_D` to *lat_mask*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing lateral boundary mask :::"
    print_text(s, cls=cls)
    self.assign_variable(self.lat_mask, lat_mask, cls=cls)
  
  def init_d_x(self, d_x, cls=None):
    r"""
    Set x-component of normalized driving-stress :math:`d_x` to *d_x*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing x-component-normalized-driving-stress direction :::"
    print_text(s, cls=cls)
    self.assign_variable(self.d_x, d_x, cls=cls)
  
  def init_d_y(self, d_y, cls=None):
    r"""
    Set y-component of normalized driving-stress :math:`d_y` to *d_y*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing y-component-normalized-driving-stress direction :::"
    print_text(s, cls=cls)
    self.assign_variable(self.d_y, d_y, cls=cls)
  
  def init_time_step(self, dt, cls=None):
    r"""
    Set time step :math:`\Delta t` to *dt*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing time step :::"
    print_text(s, cls=cls)
    self.assign_variable(self.time_step, dt, cls=cls)
  
  def init_lat(self, lat, cls=None):
    r"""
    Set grid latitude to *lat*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing grid latitude :::"
    print_text(s, cls=cls)
    self.assign_variable(self.lat, lat, cls=cls)
  
  def init_lon(self, lon, cls=None):
    r"""
    Set grid longitude to *lon*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing grid longitude :::"
    print_text(s, cls=cls)
    self.assign_variable(self.lon, lon, cls=cls)
  
  def init_M_ii(self, M_ii, cls=None):
    r"""
    Set :math:`M_{ii}` to *M_ii*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing M_ii :::"
    print_text(s, cls=cls)
    self.assign_variable(self.M_ii, M_ii, cls=cls)
  
  def init_M_ij(self, M_ij, cls=None):
    r"""
    Set :math:`M_{ij}` to *M_ij*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing M_ij :::"
    print_text(s, cls=cls)
    self.assign_variable(self.M_ij, M_ij, cls=cls)
  
  def init_M_iz(self, M_iz, cls=None):
    r"""
    Set :math:`M_{iz}` to *M_iz*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing M_iz :::"
    print_text(s, cls=cls)
    self.assign_variable(self.M_iz, M_iz, cls=cls)
  
  def init_M_ji(self, M_ji, cls=None):
    r"""
    Set :math:`M_{ji}` to *M_ji*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing M_ji :::"
    print_text(s, cls=cls)
    self.assign_variable(self.M_ji, M_ji, cls=cls)
  
  def init_M_jj(self, M_jj, cls=None):
    r"""
    Set :math:`M_{jj}` to *M_jj*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing M_jj :::"
    print_text(s, cls=cls)
    self.assign_variable(self.M_jj, M_jj, cls=cls)
  
  def init_M_jz(self, M_jz, cls=None):
    r"""
    Set :math:`M_{jz}` to *M_jz*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing M_jz :::"
    print_text(s, cls=cls)
    self.assign_variable(self.M_jz, M_jz, cls=cls)
  
  def init_M_zi(self, M_zi, cls=None):
    r"""
    Set :math:`M_{zi}` to *M_zi*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing M_zi :::"
    print_text(s, cls=cls)
    self.assign_variable(self.M_zi, M_zi, cls=cls)
  
  def init_M_zj(self, M_zj, cls=None):
    r"""
    Set :math:`M_{zj}` to *M_zj*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing M_zj :::"
    print_text(s, cls=cls)
    self.assign_variable(self.M_zj, M_zj, cls=cls)
  
  def init_M_zz(self, M_zz, cls=None):
    r"""
    Set :math:`M_{zz}` to *M_zz*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing M_zz :::"
    print_text(s, cls=cls)
    self.assign_variable(self.M_zz, M_zz, cls=cls)
  
  def init_N_ii(self, N_ii, cls=None):
    r"""
    Set membrane stress :math:`N_{ii}` to *N_ii*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing N_ii :::"
    print_text(s, cls=cls)
    self.assign_variable(self.N_ii, N_ii, cls=cls)
  
  def init_N_ij(self, N_ij, cls=None):
    r"""
    Set membrane stress :math:`N_{ij}` to *N_ij*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing N_ij :::"
    print_text(s, cls=cls)
    self.assign_variable(self.N_ij, N_ij, cls=cls)
  
  def init_N_iz(self, N_iz, cls=None):
    r"""
    Set membrane stress :math:`N_{iz}` to *N_iz*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing N_iz :::"
    print_text(s, cls=cls)
    self.assign_variable(self.N_iz, N_iz, cls=cls)
  
  def init_N_ji(self, N_ji, cls=None):
    r"""
    Set membrane stress :math:`N_{ji}` to *N_ji*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing N_ji :::"
    print_text(s, cls=cls)
    self.assign_variable(self.N_ji, N_ji, cls=cls)
  
  def init_N_jj(self, N_jj, cls=None):
    r"""
    Set membrane stress :math:`N_{jj}` to *N_jj*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing N_jj :::"
    print_text(s, cls=cls)
    self.assign_variable(self.N_jj, N_jj, cls=cls)
  
  def init_N_jz(self, N_jz, cls=None):
    r"""
    Set membrane stress :math:`N_{jz}` to *N_jz*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing N_jz :::"
    print_text(s, cls=cls)
    self.assign_variable(self.N_jz, N_jz, cls=cls)
  
  def init_N_zi(self, N_zi, cls=None):
    r"""
    Set membrane stress :math:`N_{zi}` to *N_zi*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing N_zi :::"
    print_text(s, cls=cls)
    self.assign_variable(self.N_zi, N_zi, cls=cls)
  
  def init_N_zj(self, N_zj, cls=None):
    r"""
    Set membrane stress :math:`N_{zj}` to *N_zj*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing N_zj :::"
    print_text(s, cls=cls)
    self.assign_variable(self.N_zj, N_zj, cls=cls)
  
  def init_N_zz(self, N_zz, cls=None):
    r"""
    Set membrane stress :math:`N_{zz}` to *N_zz*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing N_zz :::"
    print_text(s, cls=cls)
    self.assign_variable(self.N_zz, N_zz, cls=cls)

  def init_alpha(self, alpha, cls=None):
    r"""
    Set temperate-zone marking coefficient :math:`\alpha` to *alpha*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing temperate-zone marking coefficient :::"
    print_text(s, cls=cls)
    self.assign_variable(self.alpha, alpha, cls=cls)

  def init_alpha_int(self, alpha_int, cls=None):
    r"""
    Set vertical integral of temperate zone marking coefficient 
    :math:`\int \alpha dz` to *alpha_int*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing temperate-zone thickness :::"
    print_text(s, cls=cls)
    self.assign_variable(self.alpha_int, alpha_int, cls=cls)

  def init_Fb(self, Fb, cls=None):
    r"""
    Set basal water discharge :math:`F_b` to *Fb*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing basal-water flux :::"
    print_text(s, cls=cls)
    self.assign_variable(self.Fb, Fb, cls=cls)

  def init_Wbar(self, Wbar, cls=None):
    r"""
    Set vertically-averaged water contnet :math:`\bar{W}` to *Wbar*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing vertically averaged water content :::"
    print_text(s, cls=cls)
    self.assign_variable(self.Wbar, Wbar, cls=cls)

  def init_temp_rat(self, temp_rat, cls=None):
    r"""
    Set ratio of temperate ice :math:`\frac{1}{H}\int \alpha dz` to *temp_rat*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing temperate zone ratio :::"
    print_text(s, cls=cls)
    self.assign_variable(self.temp_rat, temp_rat, cls=cls)

  def init_Qbar(self, Qbar, cls=None):
    r"""
    Set vertically-averaged strain heat :math:`\bar{Q}` to *Qbar*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing vertically averaged strain-heat :::"
    print_text(s, cls=cls)
    self.assign_variable(self.Qbar, Qbar, cls=cls)

  def init_PE(self, PE, cls=None):
    r"""
    Set Peclet number :math:`P_e` to *PE*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing grid Peclet number :::"
    print_text(s, cls=cls)
    self.assign_variable(self.PE, PE, cls=cls)

  def init_n_f(self, n, cls=None):
    r"""
    Set outward-normal vector :math:`\mathbf{n}`` to *n*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing outward-normal-vector function n_f :::"
    print_text(s, cls=cls)
    self.assign_variable(self.n_f, n, cls=cls)

  def init_Fb_bounds(self, Fb_min, Fb_max, cls=None):
    r"""
    Set upper and lower bounds for basal water discharge :math:`F_b^{\max}` 
    and :math:`F_b^{\min}` to *Fb_max* and *Fb_min*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing bounds for basal water flux Fb :::"
    print_text(s, cls=cls)
    self.init_Fb_min(Fb_min, cls)
    self.init_Fb_max(Fb_max, cls)

  def init_Fb_min(self, Fb_min, cls=None):
    r"""
    Set lowere bound of basal water discharge :math:`F_b` to *Fb_min*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing lower bound for basal water flux Fb :::"
    print_text(s, cls=cls)
    self.assign_variable(self.Fb_min, Fb_min, cls=cls)

  def init_Fb_max(self, Fb_max, cls=None):
    r"""
    Set upper bound of basal water discharge :math:`F_b` to *Fb_max*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing upper bound for basal water flux Fb :::"
    print_text(s, cls=cls)
    self.assign_variable(self.Fb_max, Fb_max, cls=cls)

  def init_k_0(self, k_0, cls=None):
    r"""
    Set non-advective water flux coefficient :math:`k_0` to *k_0*.
    """
    if cls is None:
      cls = self.this
    s = "::: initializing non-advective flux coefficient k_0 :::"
    print_text(s, cls=cls)
    self.assign_variable(self.k_0, k_0, cls=cls)

  def init_beta_SIA(self, U_mag=None, eps=0.5):
    r"""
    Init :math:`\beta` from  :math:`\tau_b = \tau_d`, the shallow ice 
    approximation, using the observed surface velocity *U_mag* as approximate 
    basal velocity,

    .. math::

       \beta \Vert U_b \Vert = \rho g H \Vert \nabla S \Vert

    Args:

      :U_mag: basal velocity magnitude.

      :eps:   minimum velocity, default is 0.5 m/a.
    
    """
    s = "::: initializing beta from SIA :::"
    print_text(s, cls=self.this)
    Q        = self.Q
    rhoi     = self.rhoi
    g        = self.g
    gradS    = grad(self.S)
    H        = self.S - self.B
    U_s      = Function(Q, name='U_s')
    if U_mag == None:
      U_v                        = self.U_ob.vector().array()
      Ubar_v                     = self.Ubar.vector().array()
      U_v[self.Uob_missing_dofs] = Ubar_v[self.Uob_missing_dofs]
    else:
      U_v = U_mag.vector().array()
    U_v[U_v < eps] = eps
    self.assign_variable(U_s, U_v, cls=self.this)
    S_mag    = sqrt(inner(gradS, gradS) + DOLFIN_EPS)
    beta_0   = project((rhoi*g*H*S_mag) / U_s, Q, annotate=False)
    beta_0_v = beta_0.vector().array()
    beta_0_v[beta_0_v < 1e-2] = 1e-2
    self.betaSIA = Function(Q, name='betaSIA')
    self.assign_variable(self.betaSIA, beta_0_v, cls=self.this)
    
    if self.dim == 3:
      self.assign_variable(self.beta, DOLFIN_EPS, cls=self.this)
      bc_beta = DirichletBC(self.Q, self.betaSIA, self.ff, self.GAMMA_B_GND)
      bc_beta.apply(self.beta.vector())
      #self.assign_variable(self.beta, self.betaSIA, cls=self.this)
    elif self.dim == 2:
      self.assign_variable(self.beta, self.betaSIA, cls=self.this)
    print_min_max(self.beta, 'beta', cls=self.this)
      
  def init_beta_SIA_new_slide(self, U_mag=None, eps=0.5):
    r"""
    Init :math:`\beta` from  :math:`\tau_b = \tau_d`, the shallow ice 
    approximation, using the observed surface velocity *U_mag* as approximate 
    basal velocity,

    .. math::

       \beta \Vert U_b \Vert = \rho g H \Vert \nabla S \Vert

    Args:

      :U_mag: basal velocity magnitude.

      :eps:   minimum velocity, default is 0.5 m/a.
    
    """
    s = "::: initializing new sliding beta from SIA :::"
    print_text(s, cls=self.this)
    r        = 0.0
    Q        = self.Q
    rhoi     = self.rhoi
    rhow     = self.rhow
    g        = self.g
    gradS    = self.gradS
    H        = self.S - self.B
    D        = self.D
    p        = -0.383
    q        = -0.349
    
    U_s      = Function(Q, name='U_s')
    if U_mag == None:
      U_v = self.U_ob.vector().array()
    else:
      U_v = U_mag.vector().array()
    U_v[U_v < eps] = eps
    self.assign_variable(U_s, U_v, cls=self.this)
    
    Ne       = H + rhow/rhoi * D
    S_mag    = sqrt(inner(gradS, gradS) + DOLFIN_EPS)
    beta     = U_s**(1/p) / ( rhoi * g * H * S_mag * Ne**(q/p) )
    beta_0   = project(beta, Q, annotate=False)
    
    beta_0_v = beta_0.vector().array()
    beta_0_v[beta_0_v < DOLFIN_EPS] = DOLFIN_EPS
    #self.assign_variable(beta_0, beta_0_v)
    print_min_max(beta_0, 'beta_0', cls=self.this)

    #self.assign_variable(self.beta, beta_0)
    
    self.assign_variable(self.beta, DOLFIN_EPS, cls=self.this)
    bc_beta = DirichletBC(self.Q, beta_0, self.ff, GAMMA_B_GND)
    bc_beta.apply(self.beta.vector())
    
    #self.betaSIA = Function(Q)
    #self.assign_variable(self.betaSIA, beta_0_v)
    
  def init_beta_stats(self, mdl='Ubar', use_temp=False, mode='steady'):
    """
    """
    s    = "::: initializing beta from stats :::"
    print_text(s, cls=self.this)
    
    q_geo  = self.q_geo
    T_s    = self.T_surface
    adot   = self.adot
    Mb     = self.Mb
    Ubar   = self.Ubar
    Q      = self.Q
    B      = self.B
    S      = self.S
    T      = self.T
    T_s    = self.T_surface
    rho    = self.rhoi
    g      = self.g
    H      = S - B

    Ubar_v = Ubar.vector().array()
    Ubar_v[Ubar_v < 1e-10] = 1e-10
    self.assign_variable(Ubar, Ubar_v, cls=self.this)
           
    D      = Function(Q, name='D')
    B_v    = B.vector().array()
    D_v    = D.vector().array()
    D_v[B_v < 0] = B_v[B_v < 0]
    self.assign_variable(D, D_v, cls=self.this)

    gradS = as_vector([S.dx(0), S.dx(1), 0.0])
    gradB = as_vector([B.dx(0), B.dx(1), 0.0])
    gradH = as_vector([H.dx(0), H.dx(1), 0.0])

    nS   = sqrt(inner(gradS, gradS) + DOLFIN_EPS)
    nB   = sqrt(inner(gradB, gradB) + DOLFIN_EPS)
    nH   = sqrt(inner(gradH, gradH) + DOLFIN_EPS)
    
    #if mdl == 'Ubar':
    #  u_x    = -rho * g * H * S.dx(0)
    #  v_x    = -rho * g * H * S.dx(1)
    #  U_i    = as_vector([u_x,  v_x, 0.0])
    #  U_j    = as_vector([v_x, -u_x, 0.0])
    #elif mdl == 'U' or mdl == 'U_Ubar':
    #  U_i    = as_vector([self.u,  self.v, 0.0])
    #  U_j    = as_vector([self.v, -self.u, 0.0])
    U_i    = as_vector([self.u,  self.v, 0.0])
    U_j    = as_vector([self.v, -self.u, 0.0])
    Umag   = sqrt(inner(U_i,U_i) + DOLFIN_EPS)
    Uhat_i = U_i / Umag
    Uhat_j = U_j / Umag

    dBdi = dot(gradB, Uhat_i)
    dBdj = dot(gradB, Uhat_j)
    dSdi = dot(gradS, Uhat_i)
    dSdj = dot(gradS, Uhat_j)
    dHdi = dot(gradH, Uhat_i)
    dHdj = dot(gradH, Uhat_j)

    ini  = sqrt(rho * g * H * nS / (Umag + 0.1))

    x0   = S
    x1   = T_s
    x2   = nS
    x3   = D
    x4   = nB
    x5   = H
    x6   = q_geo
    x7   = adot
    x8   = T
    x9   = Mb
    x10  = self.u
    x11  = self.v
    x12  = self.w
    x13  = ln(Ubar + DOLFIN_EPS)
    x14  = ln(Umag + DOLFIN_EPS)
    x15  = ini
    x16  = dBdi
    x17  = dBdj
    x18  = dSdi
    x19  = dSdj
    x20  = nH
    x21  = self.tau_id
    x22  = self.tau_jd
    x23  = self.tau_ii
    x24  = self.tau_ij
    x25  = self.tau_ji
    x26  = self.tau_jj

    names = ['S', 'T_s', 'gradS', 'D', 'gradB', 'H', 'q_geo', 'adot', 'T',
             'Mb', 'u', 'v', 'w', 'ln(Ubar)', 'ln(Umag)', 'ini',
             'dBdi', 'dBdj', 'dSdi', 'dSdj', 'nablaH', 'tau_id', 'tau_jd',
             'tau_ii', 'tau_ij', 'tau_ji', 'tau_jj']
    names = np.array(names)

    if mdl == 'Ubar':
      if not use_temp:
        X    = [x0,x1,x5,x7,x13,x16,x18]
        idx  = [ 0, 1, 5, 7, 13, 16, 18]
        bhat = [ -1.01661102e+02,   6.59472291e-03,   8.34479667e-01,
                 -3.20751595e-04,  -1.86910058e+00,  -1.50122785e-01,
                 -1.61283407e+01,   3.42099244e+01,  -1.38190017e-07,
                 -2.42124307e-05,   5.28420031e-08,  -5.71485389e-05,
                 -3.75168897e-06,   6.62615357e-04,  -2.09616017e-03,
                 -1.63919106e-03,  -4.67468432e-07,   7.70150910e-03,
                 -1.06827565e-05,   5.82852747e-02,  -1.59176855e-01,
                  2.60703978e-08,   1.12176250e-04,  -9.96266233e-07,
                  1.54898171e-04,  -7.75201260e-03,  -3.97881378e-02,
                 -9.66212690e-04,  -6.88656946e-01,   2.86508703e+00,
                 -4.77406074e-03,   4.46234782e-03,  -9.93937326e-02,
                 -1.11058398e+01,   1.19703551e+01,  -3.46378138e+01]
        #bhat = [ -1.06707322e+02,   6.93681939e-03,   8.72090381e-01,
        #         -2.05377136e-04,  -1.68695225e+00,  -1.54427603e-01,
        #         -1.48494954e+01,   3.13320531e+01,  -1.46372911e-07,
        #         -2.54809386e-05,   5.58213888e-08,  -5.05686875e-05,
        #         -3.57485925e-06,   6.74423417e-04,  -1.90332998e-03,
        #         -1.70912922e-03,  -9.14015814e-07,   6.90894685e-03,
        #          5.38728829e-06,   5.52828014e-02,  -1.49677701e-01,
        #          2.10321794e-08,   1.26574205e-04,  -1.58804814e-06,
        #         -1.07066137e-04,  -6.59781673e-03,  -4.21221477e-02,
        #         -9.11842753e-04,  -5.91089434e-01,   2.37465616e+00,
        #         -4.79794725e-03,  -1.20787950e-03,  -8.37001425e-02,
        #         -1.35364012e+01,   2.01047113e+01,  -3.48057200e+01]
     
      else: 
        X    = [x0,x1,x5,x7,x8,x9,x13,x16,x18]
        idx  = [ 0, 1, 5, 7, 8, 9, 13, 16, 18]
        bhat = [  1.99093750e+01,  -9.37152784e-04,  -1.53849816e-03,
                 -2.72682710e-03,   3.11376629e+00,  -6.22550705e-02,
                 -4.78841821e+02,   1.18870083e-01,   1.46462501e+01,
                  4.73228083e+00,  -1.23039512e-05,   4.80948459e-08,
                 -1.75152253e-04,   1.57869882e-05,  -1.85979092e-03,
                 -5.31979350e-06,  -2.94994855e-04,  -2.88696470e-03,
                  9.87920894e-06,  -1.67014309e-02,   1.38310308e-05,
                  1.29911016e+00,   8.79462642e-06,   2.58486129e-02,
                  4.59079956e-01,  -1.62460133e-04,   8.39672735e-07,
                 -1.44977594e-02,   5.58957555e-07,   7.38625502e-04,
                 -9.92789432e-03,   6.02766800e-03,   2.74638935e-01,
                 -7.24036641e-05,  -4.63126335e-01,   2.92369712e+00,
                  5.07887934e-01,  -4.57929508e-04,  -8.33728342e-02,
                 -4.71625234e-01,  -5.85160316e-02,  -1.74723504e+01,
                 -1.83509536e+01,   5.35514345e-04,  -8.46507380e-02,
                 -1.60127263e+01]
    
    elif mdl == 'U':
      if not use_temp:
        X    = [x0,x1,x5,x7,x14,x16,x18]
        idx  = [ 0, 1, 5, 7, 14, 16, 18]
        bhat = [ -9.28289389e+01,   5.73687339e-03,   7.33526290e-01,
                  2.76998568e-03,  -1.08656857e-01,  -1.08545047e+00,
                 -1.50267782e+01,  -7.04864127e+01,  -7.76085391e-08,
                 -2.17802438e-05,  -4.99587467e-08,   5.87139196e-05,
                  1.64670170e-05,   1.06212966e-04,   7.11755177e-05,
                 -1.37677776e-03,  -9.08932836e-06,   3.60621065e-04,
                  2.97118032e-03,   5.50814766e-02,   2.21044611e-01,
                 -1.15497725e-07,   8.63993130e-05,  -2.12395318e-06,
                  7.21699958e-04,  -1.09346933e-02,  -3.12224072e-02,
                 -2.39690796e-02,  -2.95080157e-01,  -3.40502802e-01,
                 -2.62000881e-02,  -1.78157283e-02,   7.19763432e-02,
                 -1.94919730e+00,  -9.82413027e+00,  -7.61245200e+01]
      else:
        X    = [x0,x1,x5,x7,x8,x9,x14,x16,x18]
        idx  = [ 0, 1, 5, 7, 8, 9, 14, 16, 18]
        bhat = [  2.09623581e+01,   6.66919839e-04,  -7.02196170e-02,
                 -1.15080308e-03,   5.34783070e+00,  -7.11388758e-02,
                 -4.07361631e+01,   1.02018632e+00,  -1.86900651e+01,
                 -4.20181324e+01,  -9.26143019e-06,  -7.72058925e-08,
                 -4.15062408e-05,   7.02170069e-06,   2.70372865e-03,
                 -1.37333418e-05,   8.87920333e-05,   1.42938174e-03,
                  7.77557165e-06,  -2.35402146e-02,   3.04680358e-04,
                 -1.71597355e-01,   1.40252311e-04,   4.10097716e-02,
                  2.55567246e-01,  -1.33628767e-07,  -2.15459028e-06,
                  6.29599393e-05,  -4.11071912e-05,   1.28619782e-03,
                 -1.46657539e-02,   3.09279801e-03,  -2.27450062e-01,
                 -7.40025166e-03,  -5.06709113e-01,  -6.76120111e-01,
                  3.10802402e-01,  -5.34552872e-03,   2.19914707e-02,
                 -1.40943367e-01,   3.07890125e-01,  -9.03508676e+00,
                  8.27529346e+01,   6.60448755e-03,   2.42989633e+00,
                 -4.31461210e+01]
    
    elif mdl == 'U_Ubar':
      if not use_temp:
        X    = [x0,x1,x5,x7,x13,x14,x16,x18]
        idx  = [ 0, 1, 5, 7, 13, 14, 16, 18]
        bhat = [ -9.25221622e+01,   5.70295987e-03,   7.30768422e-01,
                  2.75877006e-03,   7.37861453e-02,  -2.93985236e-03,
                 -1.07390793e+00,  -1.45320123e+01,  -7.18521246e+01,
                 -7.86411913e-08,  -2.15769127e-05,  -4.80926515e-08,
                  5.56842889e-05,   1.28402687e-06,   1.12826733e-05,
                  9.07581727e-05,  -7.62357377e-05,  -1.37165484e-03,
                 -8.99331396e-06,  -3.36292037e-04,   4.24771193e-05,
                  2.97610385e-03,   5.34869351e-02,   2.28993842e-01,
                 -1.17987943e-07,   8.26468590e-05,   2.32815553e-06,
                 -6.66323072e-06,   6.73934903e-04,  -1.12192482e-02,
                 -3.22339742e-02,  -3.78492901e-04,  -2.38023512e-02,
                 -2.88687981e-01,  -4.11715791e-01,   3.06665249e-04,
                  3.29695662e-04,   4.96515338e-03,   1.28914720e-02,
                 -2.83133687e-02,  -3.08127082e-02,  -3.19074160e-02,
                 -1.60977763e+00,  -1.10451113e+01,  -7.66011531e+01]
      else:
        X    = [x0,x1,x5,x7,x8,x9,x13,x14,x16,x18]
        idx  = [ 0, 1, 5, 7, 8, 9, 13, 14, 16, 18]
        bhat = [  1.95228446e+01,   6.59477606e-04,  -6.45139002e-02,
                 -1.10071394e-03,   5.13699019e+00,  -6.45652015e-02,
                 -5.14739582e+01,  -3.68769001e-03,   9.57519905e-01,
                 -1.77507405e+01,  -4.37983921e+01,  -9.02491948e-06,
                 -7.61384926e-08,  -3.73066416e-05,   6.79516468e-06,
                  2.83564402e-03,  -4.68103812e-07,  -1.20747491e-05,
                  4.00845895e-05,   1.67755582e-03,   7.73371401e-06,
                 -2.23470170e-02,   2.78775317e-04,  -1.61211932e-01,
                  4.64633086e-05,   4.37335336e-04,   4.27466758e-02,
                  2.50573113e-01,  -4.81341231e-06,  -2.31708961e-06,
                 -1.68503900e-04,   3.54318161e-06,  -4.20165147e-05,
                  1.26878513e-03,  -1.54490818e-02,   2.66749014e-03,
                 -2.98194766e-01,  -2.92113296e-04,  -4.31378498e-03,
                 -4.83721711e-01,  -7.30055588e-01,   3.42250813e-01,
                 -3.22616161e-05,  -5.40195432e-03,   1.73408633e-02,
                 -1.31066469e-01,   9.73640123e-03,   2.61368301e-01,
                 -9.93273895e+00,   8.31773699e+01,  -5.74031885e-04,
                  9.54289863e-03,  -3.57353698e-02,   3.62295735e-03,
                  2.54399352e+00,  -4.21129483e+01]
    
    elif mdl == 'stress':
      X    = [x0,x1,x5,x7,x14,x16,x18,x21,x23,x24,x25,x26]
      idx  = [ 0, 1, 5, 7, 14, 16, 18, 21, 23, 24, 25, 26]
      bhat = [  5.47574225e+00,   9.14001489e-04,  -1.03229081e-03,
               -7.04987042e-04,   2.15686223e+00,  -1.52869679e+00,
               -1.74593819e+01,  -2.05459701e+01,  -1.23768850e-05,
                2.01460255e-05,   1.97622781e-05,   3.68067438e-05,
                6.63468606e-06,  -3.69046174e-06,  -4.47828887e-08,
               -3.67070759e-05,   2.53827543e-05,  -1.88069561e-05,
                2.05942231e-03,  -5.95566325e-10,   1.00881255e-09,
                6.11553989e-10,  -4.11737126e-10,   6.27370976e-10,
                3.42275389e-06,  -8.17017771e-03,   4.01803819e-03,
                6.78767571e-02,   4.29444354e-02,   4.45551518e-08,
               -8.23509210e-08,  -7.90182526e-08,  -1.48650850e-07,
               -2.36138203e-08,  -4.75130905e-05,  -1.81655894e-05,
                9.79852186e-04,  -1.49411705e-02,  -2.35701903e-10,
                2.32406866e-09,   1.48224703e-09,  -1.09016625e-09,
               -1.31162142e-09,   1.47593911e-02,  -1.84965301e-01,
               -1.62413731e-01,   2.38867744e-07,   2.09579112e-07,
                6.11572155e-07,   1.44891826e-06,  -4.94537953e-07,
               -3.30400642e-01,   7.93664407e-01,   7.76571489e-08,
               -1.64476914e-07,  -2.13414311e-07,   4.75810302e-07,
                2.55787543e-07,  -6.37972323e+00,  -3.77364196e-06,
                8.65062737e-08,   6.13207853e-06,   8.39233482e-07,
               -3.76402983e-06,  -2.02633500e-05,  -7.28788200e-06,
               -2.72030382e-05,  -1.33298507e-05,   1.11838930e-05,
                9.74762098e-14,  -2.37844072e-14,  -1.11310490e-13,
                8.91237008e-14,   1.16770903e-13,   5.77230478e-15,
               -4.87322338e-14,   9.62949381e-14,  -2.12122129e-13,
                1.55871983e-13]
   
    for xx,nam in zip(X, names[idx]):
      print_min_max(xx, nam, cls=self.this)

    X_i  = []
    X_i.extend(X)
     
    for i,xx in enumerate(X):
      if mdl == 'Ubar' or mdl == 'U' and not use_temp:
        k = i
      else:
        k = i+1
      for yy in X[k:]:
        X_i.append(xx*yy)
    
    #self.beta_f = exp(Constant(bhat[0]))
    self.beta_f = Constant(bhat[0])
    
    for xx,bb in zip(X_i, bhat[1:]):
      self.beta_f += Constant(bb)*xx
      #self.beta_f *= exp(Constant(bb)*xx)
    self.beta_f = exp(self.beta_f)**2
    
    #if mode == 'steady':
    #  beta0                   = project(self.beta_f, Q, annotate=False)
    #  beta0_v                 = beta0.vector().array()
    #  beta0_v[beta0_v < 1e-2] = 1e-2
    #  self.assign_variable(beta0, beta0_v, cls=self.this)
    #
    #  self.assign_variable(self.beta, 1e-2, cls=self.this)
    #  bc_beta = DirichletBC(self.Q, beta0, self.ff, self.GAMMA_B_GND)
    #  bc_beta.apply(self.beta.vector())
    
    if mode == 'steady':
      beta0  = project(self.beta_f, Q, annotate=False)
      beta0_v                 = beta0.vector().array()
      beta0_v[beta0_v < DOLFIN_EPS] = DOLFIN_EPS
      self.init_beta(beta0_v, cls=self.this)
    elif mode == 'transient':
      self.assign_variable(self.beta, 200.0, cls=self.this)
    
    print_min_max(self.beta, 'beta_hat', cls=self.this)
 
  def update_stats_beta(self):
    """
    Re-compute the statistical friction field and save into model.beta.
    """
    s    = "::: updating statistical beta :::"
    print_text(s, self.D3Model_color)
    beta   = project(self.beta_f, self.Q, annotate=False)
    beta_v = beta.vector().array()
    ##betaSIA_v = self.betaSIA.vector().array()
    ##beta_v[beta_v < 10.0]   = betaSIA_v[beta_v < 10.0]
    beta_v[beta_v < 0.0]    = 0.0
    #beta_v[beta_v > 2500.0] = 2500.0
    self.assign_variable(self.beta, beta_v, cls=self.this)
     
  def init_b_SIA(self, b, U_ob, gradS):
    r"""
    Init rate-factor b from U_ob. 
    """
    s = "::: initializing b from U_ob :::"
    print_text(s, cls=self.this)
   
    x      = self.x
    S      = self.S
    Q      = self.Q
    rhoi   = self.rhoi
    rhow   = self.rhow
    g      = self.g
    u      = U_ob[0]
    v      = U_ob[1]
    n      = 3.0
    
    class Depth(Expression):
      def eval(self, values, x):
        values[0] = min(0, x[2])
    D = Depth(element = Q.ufl_element())
    
    N      = FacetNormal(self.mesh)
    dSde   = ds(4)
           
    b_f    = TrialFunction(Q)
    phi    = TestFunction(Q)

    epi    = self.strain_rate(U_ob)
    ep_xx  = epi[0,0]
    ep_yy  = epi[1,1]
    ep_xy  = epi[0,1]
    ep_xz  = epi[0,2]
    ep_yz  = epi[1,2]   
 
    epsdot = ep_xx**2 + ep_yy**2 + ep_xx*ep_yy + ep_xy**2 + ep_xz**2 + ep_yz**2
    eta    = 0.5 * b_f * (epsdot + 1e-10)**((1-n)/(2*n))

    f_w    = rhoi*g*(S - x[2]) + rhow*g*D

    epi_1  = as_vector([   2*u.dx(0) + v.dx(1), 
                        0.5*(u.dx(1) + v.dx(0)),
                        0.5* u.dx(2)            ])
    epi_2  = as_vector([0.5*(u.dx(1) + v.dx(0)),
                             u.dx(0) + 2*v.dx(1),
                        0.5* v.dx(2)            ])

    R  = - 2 * eta * dot(epi_1, grad(phi)) * dx \
         + rhoi * g * gradS[0] * phi * dx \
         #+ 2 * eta * dot(epi_2, grad(phi)) * dx \
         #+ rhoi * g * gradS[1] * phi * dx \
   
    b_f = Function(Q)
    solve(lhs(R) == rhs(R), b_f, annotate=False)
    self.assign_variable(b, b_f, cls=self.this)

  def form_energy_dependent_rate_factor(self):
    """
    formulates energy-dependent rate factor A.
    """
    s = "::: formulating energy-dependent rate-factor :::"
    print_text(s, cls=self.this)
    
    Tp          = self.Tp
    W           = self.W
    R           = self.R
    E           = self.E
    E_shf       = self.E_shf
    E_gnd       = self.E_gnd
    a_T         = conditional( lt(Tp, 263.15),  self.a_T_l, self.a_T_u)
    Q_T         = conditional( lt(Tp, 263.15),  self.Q_T_l, self.Q_T_u)
    W_T         = conditional( lt(W,  0.01),    W,          0.01)
    self.A_shf  = E_shf*a_T*(1 + 181.25*W_T)*exp(-Q_T/(R*Tp))
    self.A_gnd  = E_gnd*a_T*(1 + 181.25*W_T)*exp(-Q_T/(R*Tp))
    self.A      = E*a_T*(1 + 181.25*W_T)*exp(-Q_T/(R*Tp))

  def calc_A(self):
    """
    calculates rate-factor A, set to model.A.
    """
    Tp          = self.Tp
    W           = self.W
    R           = self.R
    E           = self.E
    a_T         = conditional( lt(Tp, 263.15),  self.a_T_l, self.a_T_u)
    Q_T         = conditional( lt(Tp, 263.15),  self.Q_T_l, self.Q_T_u)
    W_T         = conditional( lt(W,  0.01),    W,          0.01)
    A           = E*a_T*(1 + 181.25*W_T)*exp(-Q_T/(R*Tp))
    self.A      = A
    #A           = project(A, annotate=False)
    #self.init_A(A)
 
  def calc_eta(self, epsdot):
    """
    Calculates viscosity, set to model.eta.
    """
    s     = "::: calculating viscosity :::"
    print_text(s, cls=self.this)
    eps_reg = self.eps_reg
    A       = self.A
    n       = self.n

    # calculate viscosity :
    eta     = 0.5 * A**(-1/n) * (epsdot + eps_reg)**((1-n)/(2*n))
    self.eta = eta
    #eta     = project(eta, annotate=False)
    #self.init_eta(eta)

  def calc_vert_average(self, u):
    """
    Calculates the vertical average of a given function space and function.  
    
    :param u: Function representing the model's function space
    :rtype:   Dolfin projection and Function of the vertical average
    """
    raiseNotDefined()

  def calc_normal_vector(self, cls=None):
    """
    calculates the outward normal vector as a FEniCS function.  This could then
    be used in any DirichletBC.  Saved to self.n_f.
    """
    s     = "::: calculating normal-vector function :::"
    print_text(s, cls=cls)

    n       = self.N
    n_trial = TrialFunction(self.V)
    n_test  = TestFunction(self.V)

    a = inner(n_trial, n_test)*dx
    L = inner(n,       n_test)*ds

    A = assemble(a, keep_diagonal=True)
    A.ident_zeros() # Regularize the matrix
    b = assemble(L)

    n = Function(self.V)
    solve(A, n.vector(), b, 'cg', 'amg')
    
    area = assemble(Constant(1.0)*ds(self.mesh))
    nds  = assemble(inner(n, n)*ds)
    s = "    - average value of normal on boundary: %.3f - " % (nds / area)
    print_text(s, cls=cls)
    
    self.init_n_f(n, cls=cls)

  def get_xy_velocity_angle(self, U):
    """
    Returns the angle in radians of the horizontal velocity vector from 
    the x-axis.
    """
    u,v,w   = U.split(True)
    u_v     = u.vector().array()
    v_v     = v.vector().array()
    theta_v = np.arctan2(v_v, u_v)
    Q       = u.function_space()
    theta   = Function(Q, name='theta_xy_U_angle')
    self.assign_variable(theta, theta_v, cls=self.this)
    return theta

  def get_xz_velocity_angle(self):
    """
    Returns the angle in radians of the vertical velocity vector from 
    the x-axis.
    """
    u,v,w   = self.U3.split(True)
    u_v     = u.vector().array()
    w_v     = w.vector().array()
    theta_v = np.arctan2(w_v, u_v)
    theta   = Function(self.Q, name='theta_xz_U_angle')
    self.assign_variable(theta, theta_v, cls=self.this)
    return theta

  def z_rotation_matrix(self, theta):
    """
    get rotation matrix about the z axes by angle <theta>.
    """
    c  = cos(theta)
    s  = sin(theta)
    Rz = as_matrix([[c, -s, 0],
                    [s,  c, 0],
                    [0,  0, 1]])
    return Rz

  def y_rotation_matrix(self, theta):
    """
    get rotation matrix about the y axes by angle <theta>.
    """
    c  = cos(theta)
    s  = sin(theta)
    Ry = as_matrix([[ c, 0, s],
                    [ 0, 1, 0],
                    [-s, 0, c]])
    return Ry

  def rotate_tensor(self, M, R):
    """
    rotate the tnesor <M> by the rotation matrix <R>.
    """
    if len(M.ufl_shape) == 2:
      Mr = dot(R, dot(M, R.T))
    elif len(M.ufl_shape) == 1:
      Mr = dot(R, M)
    else:
      s   = ">>> METHOD 'rotate_tensor' REQUIRES RANK 2 OR 1 TENSOR <<<"
      print_text(s, 'red', 1)
      sys.exit(1)
    return Mr

  def get_norm(self, U, type='l2'):
    """
    returns the norm of vector <U>.
    """
    # iterate through each component and convert to array :
    U_v = []
    for u in U:
      # convert to array and normailze the components of U :
      u_v = u.vector().array()
      U_v.append(u_v)
    U_v = np.array(U_v)

    # calculate the norm :
    if type == 'l2':
      norm_u = np.sqrt(np.sum(U_v**2,axis=0))
    elif type == 'linf':
      norm_u = np.amax(U_v,axis=0)
    
    return U_v, norm_u

  def normalize_vector(self, U):
    """
    Create a normalized vector of the vector <U>.
    """
    s   = "::: normalizing vector :::"
    print_text(s, cls=self.this)
    
    Q = U[0].function_space()

    U_v, norm_u = self.get_norm(U)

    norm_u[norm_u <= 0.0] = 1e-15
    
    # normalize the vector :
    U_v /= norm_u
    
    # convert back to fenics :
    U_f = []
    for u_v in U_v:
      u_f = Function(Q, name='u_f')
      self.assign_variable(u_f, u_v, cls=self.this)
      U_f.append(u_f)

    # return a UFL vector :
    return as_vector(U_f)

  def assign_submesh_variable(self, u_to, u_from):
    """
    """
    s   = "::: assigning submesh variable :::"
    print_text(s, cls=self.this)
    lg = LagrangeInterpolator()
    lg.interpolate(u_to, u_from)
    print_min_max(u_to, u_to.name(), cls=self)

  def assign_variable(self, u, var, cls=None, annotate=False):
    """
    Manually assign the values from <var> to Function <u>.  <var> may be an
    array, float, Expression, or Function.
    """
    if cls is None:
      cls = super(type(self), self)
    if isinstance(var, float) or isinstance(var, int):
      if    isinstance(u, GenericVector) or isinstance(u, Function) \
         or isinstance(u, dolfin.functions.function.Function):
        u.vector()[:] = var
      elif  isinstance(u, Constant):
        u.assign(var)
      elif  isinstance(u, float) or isinstance(u, int):
        u = var
    
    elif isinstance(var, np.ndarray):
      if var.dtype != np.float64:
        var = var.astype(np.float64)
      u.vector().set_local(var)
      u.vector().apply('insert')
    
    elif isinstance(var, Expression) \
      or isinstance(var, Constant)  \
      or isinstance(var, dolfin.functions.constant.Constant) \
      or isinstance(var, Function) \
      or isinstance(var, dolfin.functions.function.Function) \
      or isinstance(var, GenericVector):
      u.assign(var, annotate=annotate)
      #u.interpolate(var, annotate=annotate)

    #elif isinstance(var, GenericVector):
    #  self.assign_variable(u, var.array(), annotate=annotate, cls=cls)

    elif isinstance(var, str):
      File(var) >> u

    elif isinstance(var, HDF5File):
      var.read(u, u.name())

    else:
      s =  "*************************************************************\n" + \
           "assign_variable() function requires a Function, array, float,\n" + \
           " int, Vector, Expression, Constant, or string path to .xml,\n"   + \
           "not %s.  Replacing object entirely\n" + \
           "*************************************************************"
      print_text(s % type(var) , 'red', 1)
      u = var
    print_min_max(u, u.name(), cls=cls)

  def save_hdf5(self, u, f, name=None):
    """
    Save a :class:`~fenics.Function` *u* to *f* .h5 file in the hdf5 
    subdirectory of ``self.out_dir``.  If *name* = ``None``, this will 
    save the flie under ``u.name()``.
    """
    if name == None:
      name = u.name()
    s = "::: writing '%s' variable to hdf5 file :::" % name
    print_text(s, 'green')#cls=self.this)
    f.write(u, name)
    print_text("    - done -", 'green')#cls=self.this)

  def save_pvd(self, var, name, f_file=None):
    """
    Save a <name>.pvd file of the FEniCS Function <var> to this model's log 
    directory specified by model.out_dir.  If <f_file> is a File object, save 
    to this instead.
    """
    if f_file != None:
      s       = "::: saving %s.pvd file :::" % name
      print_text(s, 'green')#cls=self.this)
      f_file << var
    else:
      s       = "::: saving %spvd/%s.pvd file :::" % (self.out_dir, name)
      print_text(s, 'green')#cls=self.this)
      File(self.out_dir + 'pvd/' + name + '.pvd') << var

  def save_xdmf(self, var, name, f_file=None, t=0.0):
    """
    Save a <name>.xdmf file of the FEniCS Function <var> to this model's log 
    directory specified by model.out_dir.  If <f_file> is a File object, save 
    to this instead.
    """
    if f_file != None:
      s       = "::: saving %s.xdmf file :::" % name
      print_text(s, 'green')#cls=self.this)
      f_file << (var, float(t))
    else :
      s       = "::: saving %sxdmf/%s.xdmf file :::" % (self.out_dir, name)
      print_text(s, 'green')#cls=self.this)
      f_file = XDMFFile(self.out_dir + 'xdmf/' +  name + '.xdmf')
      f_file.write(var)
    
  def save_list_to_hdf5(self, lst, h5File):
    """
    save a list of functions or coefficients <lst> to hdf5 file <h5File>.
    """
    s    = '::: saving variables in list arg post_tmc_save_vars :::'
    print_text(s, cls=self.this)
    for var in lst:
      self.save_hdf5(var, f=h5File)

  def save_subdomain_data(self, h5File):
    """
    save all the subdomain data to hd5f file <h5File>.
    """
    raiseNotDefined()

  def save_mesh(self, h5File): 
    """
    save the mesh to hdf5 file <h5File>.
    """
    s = "::: writing 'mesh' to supplied hdf5 file :::"
    print_text(s, cls=self.this)
    h5File.write(self.mesh, 'mesh')
  
  def solve_hydrostatic_pressure(self, annotate=False, cls=None):
    """
    Solve for the hydrostatic pressure 'p'.
    """
    raiseNotDefined()
  
  def initialize_variables(self):
    """
    Initializes the class's variables to default values that are then set
    by the individually created model.
    """
    s = "::: initializing basic variables :::"
    print_text(s, cls=self.this)

    # Coordinates of various types 
    self.x             = SpatialCoordinate(self.mesh)
    self.h             = CellSize(self.mesh)
    self.N             = FacetNormal(self.mesh)
    self.lat           = Function(self.Q, name='lat')
    self.lon           = Function(self.Q, name='lon')
    self.n_f           = Function(self.V, name='n_f')

    # time step :
    self.time_step = Constant(100.0)
    self.time_step.rename('time_step', 'time step')

    # shelf mask (2 if shelf) :
    self.mask          = Function(self.Q, name='mask')
    self.init_mask(1.0) # default to all grounded ice 

    # lateral boundary mask (1 if on lateral boundary) :
    self.lat_mask      = Function(self.Q, name='lat_mask')

    # velocity mask (1 if velocity measurements present) :
    self.U_mask        = Function(self.Q, name='U_mask')

    # topography :
    self.S             = Function(self.Q_non_periodic, name='S')
    self.B             = Function(self.Q_non_periodic, name='B')
    
    # velocity observations :
    self.U_ob          = Function(self.Q, name='U_ob')
    self.u_ob          = Function(self.Q, name='u_ob')
    self.v_ob          = Function(self.Q, name='v_ob')
    
    # unified velocity (non-periodic because it is always known everywhere) :
    self.U_mag         = Function(self.Q,               name='U_mag')
    self.U3            = Function(self.Q3_non_periodic, name='U3')
    u,v,w              = self.U3.split()
    u.rename('u', '')
    v.rename('v', '')
    w.rename('w', '')
    self.u             = u
    self.v             = v
    self.w             = w
    self.assx  = FunctionAssigner(u.function_space(), self.Q_non_periodic)
    self.assy  = FunctionAssigner(v.function_space(), self.Q_non_periodic)
    self.assz  = FunctionAssigner(w.function_space(), self.Q_non_periodic)

    # momentum model :
    self.eta           = Function(self.Q, name='eta')
    self.p             = Function(self.Q_non_periodic, name='p')
    self.beta          = Function(self.Q, name='beta')
    self.E             = Function(self.Q, name='E')
    self.E_gnd         = Function(self.Q, name='E_gnd')
    self.E_shf         = Function(self.Q, name='E_shf')
    self.A             = Function(self.Q, name='A')
    self.A_gnd         = Function(self.Q, name='A_gnd')
    self.A_shf         = Function(self.Q, name='A_shf')
    self.u_lat         = Function(self.Q, name='u_lat')
    self.v_lat         = Function(self.Q, name='v_lat')
    self.w_lat         = Function(self.Q, name='w_lat')
    self.lam           = Function(self.Q, name='lam')
    self.init_E(1.0) # always use no enhancement on rate-factor A 
    
    # energy model :
    self.T             = Function(self.Q, name='T')
    self.Tp            = Function(self.Q, name='Tp')
    self.q_geo         = Function(self.Q, name='q_geo')
    self.q_fric        = Function(self.Q, name='q_fric')
    self.gradT_B       = Function(self.Q, name='gradT_B')
    self.gradTm_B      = Function(self.Q, name='gradTm_B')
    self.theta         = Function(self.Q, name='theta')
    self.W             = Function(self.Q, name='W')
    self.Wc            = Function(self.Q, name='Wc')
    self.Mb            = Function(self.Q, name='Mb')
    self.rhob          = Function(self.Q, name='rhob')
    self.T_melt        = Function(self.Q, name='T_melt')     # pressure-melting
    self.theta_melt    = Function(self.Q, name='theta_melt') # pressure-melting
    self.T_surface     = Function(self.Q, name='T_surface')
    self.alpha         = Function(self.Q, name='alpha')
    self.alpha_int     = Function(self.Q, name='alpha_int')
    self.Fb            = Function(self.Q, name='Fb')
    self.PE            = Function(self.Q, name='PE')
    self.Wbar          = Function(self.Q, name='Wbar')
    self.Fb_min        = Function(self.Q, name='Fb_min')
    self.Fb_max        = Function(self.Q, name='Fb_max')
    self.Qbar          = Function(self.Q, name='Qbar')
    self.temp_rat      = Function(self.Q, name='temp_rat')
    self.k_0           = Constant(1.0,    name='k_0')
    self.k_0.rename('k_0', 'k_0')
    
    # adjoint model :
    self.control_opt   = Function(self.Q, name='control_opt')

    # balance Velocity model :
    self.adot          = Function(self.Q, name='adot')
    self.d_x           = Function(self.Q, name='d_x')
    self.d_y           = Function(self.Q, name='d_y')
    self.Ubar          = Function(self.Q, name='Ubar')
    self.uhat          = Function(self.Q, name='uhat')
    self.vhat          = Function(self.Q, name='vhat')
    
    # Stress-balance model (this is always non-periodic) :
    self.M_ii          = Function(self.Q_non_periodic, name='M_ii')
    self.M_ij          = Function(self.Q_non_periodic, name='M_ij')
    self.M_iz          = Function(self.Q_non_periodic, name='M_iz')
    self.M_ji          = Function(self.Q_non_periodic, name='M_ji')
    self.M_jj          = Function(self.Q_non_periodic, name='M_jj')
    self.M_jz          = Function(self.Q_non_periodic, name='M_jz')
    self.M_zi          = Function(self.Q_non_periodic, name='M_zi')
    self.M_zj          = Function(self.Q_non_periodic, name='M_zj')
    self.M_zz          = Function(self.Q_non_periodic, name='M_zz')
    self.N_ii          = Function(self.Q_non_periodic, name='N_ii')
    self.N_ij          = Function(self.Q_non_periodic, name='N_ij')
    self.N_iz          = Function(self.Q_non_periodic, name='N_iz')
    self.N_ji          = Function(self.Q_non_periodic, name='N_ji')
    self.N_jj          = Function(self.Q_non_periodic, name='N_jj')
    self.N_jz          = Function(self.Q_non_periodic, name='N_jz')
    self.N_zi          = Function(self.Q_non_periodic, name='N_zi')
    self.N_zj          = Function(self.Q_non_periodic, name='N_zj')
    self.N_zz          = Function(self.Q_non_periodic, name='N_zz')

  def home_rolled_newton_method(self, R, U, J, bcs, atol=1e-7, rtol=1e-10,
                                relaxation_param=1.0, max_iter=25,
                                method='cg', preconditioner='amg',
                                cb_ftn=None):
    """
    Appy Newton's method.

    Args:
     
      :R:                residual of system
      :U:                unknown to determine
      :J:                Jacobian
      :bcs:              set of Dirichlet boundary conditions
      :atol:             absolute stopping tolerance
      :rtol:             relative stopping tolerance
      :relaxation_param: ratio of down-gradient step to take each iteration.
      :max_iter:         maximum number of iterations to perform
      :method:           linear solution method
      :preconditioner:   preconditioning method to use with ``Krylov`` solver
      :cb_ftn:           at the end of each iteration, this is called
    """
    converged  = False
    lmbda      = relaxation_param   # relaxation parameter
    nIter      = 0                  # number of iterations
   
    # need to homogenize the boundary, as the residual is always zero over
    # essential boundaries :
    bcs_u = []
    for bc in bcs:
      bc = DirichletBC(bc)
      bc.homogenize()
      bcs_u.append(bc)
    
    # the direction of decent :
    d = Function(U.function_space()) 
    
    while not converged and nIter < max_iter:
      # assemble system :
      A, b    = assemble_system(J, -R, bcs_u)
      
      # determine step direction :
      solve(A, d.vector(), b, method, preconditioner, annotate=False)
    
      # calculate residual :
      residual  = b.norm('l2')
    
      # set initial residual : 
      if nIter == 0:
        residual_0 = residual

      # the relative residual :
      rel_res = residual/residual_0

      # check for convergence :
      converged = residual < atol or rel_res < rtol
      
      # move U down the gradient :
      U.vector()[:] += lmbda*d.vector()
      
      # increment counter :
      nIter += 1
    
      # print info to screen :
      if self.MPI_rank == 0:
        string = "Newton iteration %d: r (abs) = %.3e (tol = %.3e) " \
                 +"r (rel) = %.3e (tol = %.3e)"
        print string % (nIter, residual, atol, rel_res, rtol)

      # call the callback function, if desired :
      if cb_ftn is not None:
        s    = "::: calling home-rolled Newton method callback :::"
        print_text(s, cls=self.this)
        cb_ftn()

  def thermo_solve(self, momentum, energy, wop_kwargs,
                   callback=None, atol=1e2, rtol=1e0, max_iter=50,
                   iter_save_vars=None, post_tmc_save_vars=None,
                   starting_i=1):
    r""" 
    Perform thermo-mechanical coupling between momentum and energy.

    Args:

      :momentum:       a :class:`~momentum.Momentum` instance
      :energy:         a :class:`~energy.Energy` instance.  Currently this only 
                       works for :class:`~energy.Enthalpy`
      :wop_kwargs:     a :py:class:`~dict` of arguments for water-optimization 
                       method :func:`~energy.Energy.optimize_water_flux`
      :callback:       a function that is called back at the end of each 
                       iteration
      :atol:           absolute stopping tolerance 
                       :math:`a_{tol} \leq r = \Vert \theta_n - \theta_{n-1} \Vert`
      :rtol:           relative stopping tolerance
                       :math:`r_{tol} \leq \Vert r_n - r_{n-1} \Vert`
      :max_iter:       maximum number of iterations to perform
      :iter_save_vars: python :py:class:`~list` containing functions to 
                       save each iteration
      :starting_i:     if you are restarting this process, you may start 
                       it at a later iteration. 
    """
    s    = '::: performing thermo-mechanical coupling with atol = %.2e, ' + \
           'rtol = %.2e, and max_iter = %i :::'
    print_text(s % (atol, rtol, max_iter), cls=self.this)
    
    from cslvr import Momentum
    from cslvr import Energy
    
    if not isinstance(momentum, Momentum):
      s = ">>> thermo_solve REQUIRES A 'Momentum' INSTANCE, NOT %s <<<"
      print_text(s % type(momentum) , 'red', 1)
      sys.exit(1)
    
    if not isinstance(energy, Energy):
      s = ">>> thermo_solve REQUIRES AN 'Energy' INSTANCE, NOT %s <<<"
      print_text(s % type(energy) , 'red', 1)
      sys.exit(1)

    # mark starting time :
    t0   = time()

    # ensure that we have a steady-state form :
    if energy.transient:
      energy.make_steady_state()

    # retain base install directory :
    out_dir_i = self.out_dir

    # directory for saving convergence history :
    d_hist   = self.out_dir + 'tmc/convergence_history/'
    if not os.path.exists(d_hist) and self.MPI_rank == 0:
      os.makedirs(d_hist)

    # number of digits for saving variables :
    n_i  = len(str(max_iter))
    
    # get the bounds of Fb, the max will be updated based on temperate zones :
    if energy.energy_flux_mode == 'Fb':
      bounds = copy(wop_kwargs['bounds'])
      self.init_Fb_min(bounds[0], cls=self.this)
      self.init_Fb_max(bounds[1], cls=self.this)
      wop_kwargs['bounds']  = (self.Fb_min, self.Fb_max)

    # L_2 erro norm between iterations :
    abs_error = np.inf
    rel_error = np.inf
      
    # number of iterations, from a starting point (useful for restarts) :
    if starting_i <= 1:
      counter = 1
    else:
      counter = starting_i
   
    # previous velocity for norm calculation
    U_prev    = self.theta.copy(True)

    # perform a fixed-point iteration until the L_2 norm of error 
    # is less than tolerance :
    while abs_error > atol and rel_error > rtol and counter <= max_iter:
       
      # set a new unique output directory :
      out_dir_n = 'tmc/%0*d/' % (n_i, counter)
      self.set_out_dir(out_dir_i + out_dir_n)
      
      # solve velocity :
      momentum.solve(annotate=False)

      # update pressure-melting point :
      energy.calc_T_melt(annotate=False)

      # calculate basal friction heat flux :
      momentum.calc_q_fric()
      
      # derive temperature and temperature-melting flux terms :
      energy.calc_basal_temperature_flux()
      energy.calc_basal_temperature_melting_flux()

      # solve energy steady-state equations to derive temperate zone :
      energy.derive_temperate_zone(annotate=False)
      
      # fixed-point interation for thermal parameters and discontinuous 
      # properties :
      energy.update_thermal_parameters(annotate=False)
      
      # calculate the basal-melting rate :
      energy.solve_basal_melt_rate()
      
      # always initialize Fb to the zero-energy-flux bc :  
      Fb_v = self.Mb.vector().array() * self.rhoi(0) / self.rhow(0)
      self.init_Fb(Fb_v)
  
      # update bounds based on temperate zone :
      if energy.energy_flux_mode == 'Fb':
        Fb_m_v                 = self.Fb_max.vector().array()
        alpha_v                = self.alpha.vector().array()
        Fb_m_v[:]              = DOLFIN_EPS
        Fb_m_v[alpha_v == 1.0] = bounds[1]
        self.init_Fb_max(Fb_m_v, cls=self.this)
      
      # optimize the flux of water to remove abnormally high water :
      if energy.energy_flux_mode == 'Fb':
        energy.optimize_water_flux(**wop_kwargs)

      # solve the energy-balance and partition T and W from theta :
      energy.solve(annotate=False)
      
      # calculate L_2 norms :
      abs_error_n  = norm(U_prev.vector() - self.theta.vector(), 'l2')
      tht_nrm      = norm(self.theta.vector(), 'l2')

      # save convergence history :
      if counter == 1:
        rel_error  = abs_error_n
        if self.MPI_rank == 0:
          err_a = np.array([abs_error_n])
          nrm_a = np.array([tht_nrm])
          np.savetxt(d_hist + 'abs_err.txt',    err_a)
          np.savetxt(d_hist + 'theta_norm.txt', nrm_a)
      else:
        rel_error = abs(abs_error - abs_error_n)
        if self.MPI_rank == 0:
          err_n = np.loadtxt(d_hist + 'abs_err.txt')
          nrm_n = np.loadtxt(d_hist + 'theta_norm.txt')
          err_a = np.append(err_n, np.array([abs_error_n]))
          nrm_a = np.append(nrm_n, np.array([tht_nrm]))
          np.savetxt(d_hist + 'abs_err.txt',     err_a)
          np.savetxt(d_hist + 'theta_norm.txt',  nrm_a)

      # print info to screen :
      if self.MPI_rank == 0:
        s0    = '>>> '
        s1    = 'TMC fixed-point iteration %i (max %i) done: ' \
                 % (counter, max_iter)
        s2    = 'r (abs) = %.2e ' % abs_error
        s3    = '(tol %.2e), '    % atol
        s4    = 'r (rel) = %.2e ' % rel_error
        s5    = '(tol %.2e)'      % rtol
        s6    = ' <<<'
        text0 = get_text(s0, 'red', 1)
        text1 = get_text(s1, 'red')
        text2 = get_text(s2, 'red', 1)
        text3 = get_text(s3, 'red')
        text4 = get_text(s4, 'red', 1)
        text5 = get_text(s5, 'red')
        text6 = get_text(s6, 'red', 1)
        print text0 + text1 + text2 + text3 + text4 + text5 + text6
      
      # update error stuff and increment iteration counter :
      abs_error    = abs_error_n
      U_prev       = self.theta.copy(True)
      counter     += 1

      # call callback function if set :
      if callback != None:
        s    = '::: calling thermo-couple-callback function :::'
        print_text(s, cls=self.this)
        callback()
    
      # save state to unique hdf5 file :
      if isinstance(iter_save_vars, list):
        s    = '::: saving variables in list arg iter_save_vars :::'
        print_text(s, cls=self.this)
        out_file = self.out_dir + 'tmc.h5'
        foutput  = HDF5File(mpi_comm_world(), out_file, 'w')
        for var in iter_save_vars:
          self.save_hdf5(var, f=foutput)
        foutput.close()
    
    # reset the base directory ! :
    self.set_out_dir(out_dir_i)
    
    # reset the bounds on Fb :
    if energy.energy_flux_mode == 'Fb':  wop_kwargs['bounds'] = bounds
      
    # save state to unique hdf5 file :
    if isinstance(post_tmc_save_vars, list):
      s    = '::: saving variables in list arg post_tmc_save_vars :::'
      print_text(s, cls=self.this)
      out_file = self.out_dir + 'tmc.h5'
      foutput  = HDF5File(mpi_comm_world(), out_file, 'w')
      for var in post_tmc_save_vars:
        self.save_hdf5(var, f=foutput)
      foutput.close()

    # calculate total time to compute
    tf = time()
    s  = tf - t0
    m  = s / 60.0
    h  = m / 60.0
    s  = s % 60
    m  = m % 60
    text = "time to thermo-couple: %02d:%02d:%02d" % (h,m,s)
    print_text(text, 'red', 1)
       
    # plot the convergence history : 
    s    = "::: convergence info saved to \'%s\' :::"
    print_text(s % d_hist, cls=self.this)
    if self.MPI_rank == 0:
      np.savetxt(d_hist + 'time.txt', np.array([tf - t0]))

      err_a = np.loadtxt(d_hist + 'abs_err.txt')
      nrm_a = np.loadtxt(d_hist + 'theta_norm.txt')
     
      # plot iteration error : 
      fig   = plt.figure()
      ax    = fig.add_subplot(111)
      ax.set_ylabel(r'$\Vert \theta_{n-1} - \theta_n \Vert$')
      ax.set_xlabel(r'iteration')
      ax.plot(err_a, 'k-', lw=2.0)
      plt.grid()
      plt.savefig(d_hist + 'abs_err.png', dpi=100)
      plt.close(fig)
      
      # plot theta norm :
      fig = plt.figure()
      ax  = fig.add_subplot(111)
      ax.set_ylabel(r'$\Vert \theta_n \Vert$')
      ax.set_xlabel(r'iteration')
      ax.plot(nrm_a, 'k-', lw=2.0)
      plt.grid()
      plt.savefig(d_hist + 'theta_norm.png', dpi=100)
      plt.close(fig)

  def assimilate_U_ob(self, momentum, beta_i, max_iter, 
                      tmc_kwargs, uop_kwargs,
                      atol                = 1e2,
                      rtol                = 1e0, 
                      initialize          = True,
                      incomplete          = True,
                      post_iter_save_vars = None,
                      post_ini_callback   = None,
                      starting_i          = 1):
    """
    """
    s    = '::: performing assimilation process with %i max iterations :::'
    print_text(s % max_iter, cls=self.this)

    # retain base install directory :
    out_dir_i = self.out_dir
    
    # directory for saving convergence history :
    d_hist   = self.out_dir + 'convergence_history/'
    if not os.path.exists(d_hist) and self.MPI_rank == 0:
      os.makedirs(d_hist)

    # number of digits for saving variables :
    n_i  = len(str(max_iter))
    
    # starting time :
    t0   = time()
    
    # L_2 erro norm between iterations :
    abs_error = np.inf
    rel_error = np.inf
      
    # number of iterations, from a starting point (useful for restarts) :
    if starting_i <= 1:
      counter = 1
    else:
      counter = starting_i

    # initialize friction field :
    self.init_beta(beta_i, cls=self.this)
   
    # previous friction for norm calculation :
    beta_prev    = self.beta.copy(True)

    # perform initialization step if desired :
    if initialize:
      s    = '    - performing initialization step -'
      print_text(s, cls=self.this)

      # set the initialization output directory :
      out_dir_n = 'initialization/'
      self.set_out_dir(out_dir_i + out_dir_n)
      
      # thermo-mechanical couple :
      self.thermo_solve(**tmc_kwargs)

      # call the post function if set :
      if post_ini_callback is not None:
        s    = '::: calling post-initialization assimilate_U_ob ' + \
               'callback function :::'
        print_text(s, cls=self.this)
        post_ini_callback()

    # otherwise, tell us that we are not initializing :
    else:
      s    = '    - skipping initialization step -'
      print_text(s, cls=self.this)
    
    # save the w_opt bounds on Fb :
    bounds = copy(tmc_kwargs['wop_kwargs']['bounds'])

    # assimilate the data : 
    while abs_error > atol and rel_error > rtol and counter <= max_iter:
      s    = '::: entering iterate %i of %i of assimilation process :::'
      print_text(s % (counter, max_iter), cls=self.this)
       
      # set a new unique output directory :
      out_dir_n = '%0*d/' % (n_i, counter)
      self.set_out_dir(out_dir_i + out_dir_n)
   
      # the incomplete adjoint means the viscosity is linear, and
      # we do not want to reset the original momentum configuration, because
      # we have more non-linear solves to do :
      if incomplete and not momentum.linear:
        momentum.linearize_viscosity(reset_orig_config=True)
    
      # re-initialize friction field :
      if counter > starting_i: self.init_beta(beta_i, cls=self.this)

      # optimize the velocity : 
      momentum.optimize_U_ob(**uop_kwargs)

      # reset the momentum to the original configuration : 
      if not momentum.linear_s and momentum.linear: momentum.reset()

      # thermo-mechanically couple :
      self.thermo_solve(**tmc_kwargs)
      
      # calculate L_2 norms :
      abs_error_n  = norm(beta_prev.vector() - self.beta.vector(), 'l2')
      beta_nrm     = norm(self.beta.vector(), 'l2')

      # save convergence history :
      if counter == 1:
        rel_error  = abs_error_n
        if self.MPI_rank == 0:
          err_a = np.array([abs_error_n])
          nrm_a = np.array([beta_nrm])
          np.savetxt(d_hist + 'abs_err.txt',   err_a)
          np.savetxt(d_hist + 'beta_norm.txt', nrm_a)
      else:
        rel_error = abs(abs_error - abs_error_n)
        if self.MPI_rank == 0:
          err_n = np.loadtxt(d_hist + 'abs_err.txt')
          nrm_n = np.loadtxt(d_hist + 'beta_norm.txt')
          err_a = np.append(err_n, np.array([abs_error_n]))
          nrm_a = np.append(nrm_n, np.array([beta_nrm]))
          np.savetxt(d_hist + 'abs_err.txt',    err_a)
          np.savetxt(d_hist + 'beta_norm.txt',  nrm_a)

      # print info to screen :
      if self.MPI_rank == 0:
        s0    = '>>> '
        s1    = 'U_ob assimilation iteration %i (max %i) done: ' \
                % (counter, max_iter)
        s2    = 'r (abs) = %.2e ' % abs_error
        s3    = '(tol %.2e), '    % atol
        s4    = 'r (rel) = %.2e ' % rel_error
        s5    = '(tol %.2e)'      % rtol
        s6    = ' <<<'
        text0 = get_text(s0, 'red', 1)
        text1 = get_text(s1, 'red')
        text2 = get_text(s2, 'red', 1)
        text3 = get_text(s3, 'red')
        text4 = get_text(s4, 'red', 1)
        text5 = get_text(s5, 'red')
        text6 = get_text(s6, 'red', 1)
        print text0 + text1 + text2 + text3 + text4 + text5 + text6
      
      # save state to unique hdf5 file :
      if isinstance(post_iter_save_vars, list):
        s    = '::: saving variables in list arg post_iter_save_vars :::'
        print_text(s, cls=self.this)
        out_file = self.out_dir + 'inverted.h5'
        foutput  = HDF5File(mpi_comm_world(), out_file, 'w')
        for var in post_iter_save_vars:
          self.save_hdf5(var, f=foutput)
        foutput.close()
      
      # update error stuff and increment iteration counter :
      abs_error    = abs_error_n
      beta_prev    = self.beta.copy(True)
      counter     += 1

    # calculate total time to compute
    tf = time()
    s  = tf - t0
    m  = s / 60.0
    h  = m / 60.0
    s  = s % 60
    m  = m % 60
    text = "time to compute TMC optimized ||u - u_ob||: %02d:%02d:%02d"
    print_text(text % (h,m,s) , 'red', 1)
       
    # plot the convergence history : 
    s    = "::: convergence info saved to \'%s\' :::"
    print_text(s % d_hist, cls=self.this)
    if self.MPI_rank == 0:
      np.savetxt(d_hist + 'time.txt', np.array([tf - t0]))

      err_a = np.loadtxt(d_hist + 'abs_err.txt')
      nrm_a = np.loadtxt(d_hist + 'beta_norm.txt')
     
      # plot iteration error : 
      fig   = plt.figure()
      ax    = fig.add_subplot(111)
      ax.set_ylabel(r'$\Vert \beta_{n-1} - \beta_n \Vert$')
      ax.set_xlabel(r'iteration')
      ax.plot(err_a, 'k-', lw=2.0)
      plt.grid()
      plt.savefig(d_hist + 'abs_err.png', dpi=100)
      plt.close(fig)
      
      # plot theta norm :
      fig = plt.figure()
      ax  = fig.add_subplot(111)
      ax.set_ylabel(r'$\Vert \beta_n \Vert$')
      ax.set_xlabel(r'iteration')
      ax.plot(nrm_a, 'k-', lw=2.0)
      plt.grid()
      plt.savefig(d_hist + 'beta_norm.png', dpi=100)
      plt.close(fig)

  def L_curve(self, alphas, physics, control, int_domain, adj_ftn, adj_kwargs,
              reg_kind='Tikhonov', pre_callback=None, post_callback=None,
              itr_save_vars=None):
    """
    """
    s    = '::: starting L-curve procedure :::'
    print_text(s, cls=self.this)
    
    # starting time :
    t0   = time()

    # retain base install directory :
    out_dir_i = self.out_dir

    # retain initial control parameter for consistency :
    control_ini = control.copy(True)

    # functional lists to be populated :
    Js     = []
    Rs     = []
   
    # iterate through each of the regularization parameters provided : 
    for i,alpha in enumerate(alphas):
      s    = '::: performing L-curve iteration %i with alpha = %.3e :::'
      print_text(s % (i,alpha) , atrb=1, cls=self.this)

      # reset everything after the first iteration :
      if i > 0:
        s    = '::: initializing physics :::'
        print_text(s, cls=self.this)
        physics.reset()
        self.assign_variable(control, control_ini, cls=self.this)
      
      # set the appropriate output directory :
      out_dir_n = 'alpha_%.1E/' % alpha
      self.set_out_dir(out_dir_i + out_dir_n)
      
      # call the pre-adjoint callback function :
      if pre_callback is not None:
        s    = '::: calling L_curve() pre-adjoint pre_callback() :::'
        print_text(s, cls=self.this)
        pre_callback()
     
      # get new regularization functional : 
      R = physics.form_reg_ftn(control, integral=int_domain,
                               kind=reg_kind, alpha=alpha)

      # solve the adjoint system :
      adj_ftn(**adj_kwargs)
      
      # calculate functionals of interest :
      Rs.append(assemble(physics.Rp))
      Js.append(assemble(physics.Jp))
      
      # call the pre-adjoint callback function :
      if post_callback is not None:
        s    = '::: calling L_curve() post-adjoint post_callback() :::'
        print_text(s, cls=self.this)
        post_callback()
      
      # save state to unique hdf5 file :
      if isinstance(itr_save_vars, list):
        s    = '::: saving variables in list arg itr_save_vars :::'
        print_text(s, cls=self.this)
        out_file = self.out_dir + 'lcurve.h5'
        foutput  = HDF5File(mpi_comm_world(), out_file, 'w')
        for var in itr_save_vars:
          self.save_hdf5(var, f=foutput)
        foutput.close()
    
    s    = '::: L-curve procedure complete :::'
    print_text(s, cls=self.this)

    # calculate total time to compute
    s = time() - t0
    m = s / 60.0
    h = m / 60.0
    s = s % 60
    m = m % 60
    text = "time to complete L-curve procedure: %02d:%02d:%02d" % (h,m,s)
    print_text(text, 'red', 1)
   
    # save the resulting functional values and alphas to CSF : 
    if self.MPI_rank==0:
      # iterate through the directiories we just created and grab the data :
      alphas = []
      Ds     = []
      Js     = []
      J1s    = []
      J2s    = []
      Rs     = []
      ns     = [] 
      for d in next(os.walk(out_dir_i))[1]:
        m = re.search('(alpha_)(\d\W\dE\W\d+)', d)
        if m is not None:
          do = out_dir_i + d + '/objective_ftnls_history/'
          alphas.append(float(m.group(2)))
          Ds.append(np.loadtxt(do  + 'Ds.txt'))
          Js.append(np.loadtxt(do  + 'Js.txt'))
          J1s.append(np.loadtxt(do + 'J1s.txt'))
          J2s.append(np.loadtxt(do + 'J2s.txt'))
          Rs.append(np.loadtxt(do  + 'Rs.txt'))
          ns.append(len(Js[-1]))
      alphas = np.array(alphas) 
      Ds     = np.array(Ds) 
      Js     = np.array(Js) 
      J1s    = np.array(J1s) 
      J2s    = np.array(J2s) 
      Rs     = np.array(Rs) 
      ns     = np.array(ns)

      # sort everything :
      idx    = np.argsort(alphas)
      alphas = alphas[idx]
      Ds     = Ds[idx]
      Js     = Js[idx]
      J1s    = J1s[idx]
      J2s    = J2s[idx]
      Rs     = Rs[idx]
      ns     = ns[idx]
     
      # plot the functionals : 
      #=========================================================================
      fig = plt.figure(figsize=(6,2.5))
      ax  = fig.add_subplot(111)

      # we want to plot the different alpha values a different shade :
      cmap = plt.get_cmap('viridis')
      colors = [ cmap(x) for x in np.linspace(0, 1, len(alphas)) ]
      
      k    = 0    # counter so we can plot side-by-side
      ints = [0]  # to modify the x-axis labels
      for i,c in zip(range(len(alphas)), colors):
        xi = np.arange(k, k + ns[i])
        ints.append(xi.max())
        # if this is the first iteration, we put a legend on it :
        if i == 0:
          # if we have two cost functionals, let's plot both :
          if physics.obj_ftn_type == 'log_L2_hybrid':
            ax.plot(xi, J1s[i], '-',  c='0.5', lw=2.0,
                    label = r'$\mathscr{I}_1$')
            ax.plot(xi, J2s[i], '-',  c='k',   lw=2.0,
                    label = r'$\mathscr{I}_2$')
          # otherwise, just the one :
          else:
            ax.plot(xi, Js[i], '-',   c='k',   lw=2.0,
                    label = r'$\mathscr{I}$')
          # always plot the regularization functional :
          ax.plot(xi, Rs[i],   '-',   c='r',   lw=2.0,
                  label = r'$\mathscr{R}$')
        # otherwise, we don't need cluttered legends :
        else:
          # if we have two cost functionals, let's plot both :
          if physics.obj_ftn_type == 'log_L2_hybrid':
            ax.plot(xi, J1s[i], '-',  c='0.5', lw=2.0)
            ax.plot(xi, J2s[i], '-',  c='k',   lw=2.0)
          # otherwise, just the one :
          else:  
            ax.plot(xi, Js[i], '-',  c='k', lw=2.0)
          # always plot the regularization functional :
          ax.plot(xi, Rs[i],   '-',  c='r', lw=2.0)
        k += ns[i] - 1
      ints = np.array(ints)
      
      label = []
      for i in alphas:
        label.append(r'$\gamma = %g$' % i)
      
      # reset the x-label to be meaningfull :
      ax.set_xticks(ints)
      ax.set_xticklabels(label, size='small', ha='left')#, rotation=-45)
      ax.set_xlabel(r'relative iteration')
      
      ax.grid()
      ax.set_yscale('log')
      
      # plot the functional legend across the top in a row : 
      if physics.obj_ftn_type == 'log_L2_hybrid': ncol = 3
      else :                                      ncol = 2 
      leg = ax.legend(loc='upper center', ncol=ncol)
      leg.get_frame().set_alpha(0.0)
      
      plt.tight_layout()
      plt.savefig(out_dir_i + 'convergence.pdf')
      plt.close(fig)

      # plot L-curve :
      #=========================================================================

      # we only want the last value of each optimization : 
      fin_Js = Js[:,-1]
      fin_Rs = Rs[:,-1]
      
      fig = plt.figure(figsize=(6,2.5))
      ax  = fig.add_subplot(111)
      
      ax.plot(fin_Js, fin_Rs, 'k-', lw=2.0)
      
      ax.grid()
    
      # useful for figuring out what reg. parameter goes with what :  
      for i,c in zip(range(len(alphas)), colors):
        ax.plot(fin_Js[i], fin_Rs[i], 'o',  c=c, lw=2.0,
                label = r'$\gamma = %g$' % alphas[i])
     
      ax.set_xlabel(r'$\mathscr{I}^*$')
      ax.set_ylabel(r'$\mathscr{R}^*$')
      
      leg = ax.legend(loc='upper right', ncol=2)
      leg.get_frame().set_alpha(0.0)
      
      ax.set_yscale('log')
      #ax.set_xscale('log')
      
      plt.tight_layout()
      plt.savefig(out_dir_i + 'l_curve.pdf')
      plt.close(fig)

      # save the functionals :
      #=========================================================================

      d = out_dir_i + 'functionals/'
      if not os.path.exists(d):
        os.makedirs(d)
      np.savetxt(d + 'Rs.txt',   np.array(fin_Rs))
      np.savetxt(d + 'Js.txt',   np.array(fin_Js))
      np.savetxt(d + 'as.txt',   np.array(alphas))

  def transient_solve(self, momentum, energy, mass, t_start, t_end, time_step,
                      adaptive=False, annotate=False, callback=None):
    """
    """
    s    = '::: performing transient run :::'
    print_text(s, cls=self)
    
    from cslvr.momentum import Momentum
    from cslvr.energy   import Energy
    from cslvr.mass     import Mass
    
    if momentum.__class__.__base__ != Momentum:
      s = ">>> transient_solve REQUIRES A 'Momentum' INSTANCE, NOT %s <<<"
      print_text(s % type(momentum), 'red', 1)
      sys.exit(1)
    
    if energy.__class__.__base__ != Energy:
      s = ">>> transient_solve REQUIRES AN 'Energy' INSTANCE, NOT %s <<<"
      print_text(s % type(energy), 'red', 1)
      sys.exit(1)
    
    if mass.__class__.__base__ != Mass:
      s = ">>> transient_solve REQUIRES A 'Mass' INSTANCE, NOT %s <<<"
      print_text(s % type(mass), 'red', 1)
      sys.exit(1)
    
    stars = "*****************************************************************"
    self.init_time_step(time_step)
    self.step_time = []
    t0             = time()
    t              = t_start
    dt             = time_step
    alpha          = momentum.solve_params['solver']['newton_solver']
    alpha          = alpha['relaxation_parameter']
   
    # Loop over all times
    while t <= t_end:

      # start the timer :
      tic = time()
      
      # solve momentum equation, lower alpha on failure :
      if adaptive:
        solved_u = False
        par    = momentum.solve_params['solver']['newton_solver']
        while not solved_u:
          if par['relaxation_parameter'] < 0.2:
            status_u = [False, False]
            break
          # always reset velocity for good convergence :
          self.assign_variable(momentum.get_U(), DOLFIN_EPS, cls=self)
          status_u = momentum.solve(annotate=annotate)
          solved_u = status_u[1]
          if not solved_u:
            par['relaxation_parameter'] /= 1.43
            print_text(stars, 'red', 1)
            s = ">>> WARNING: newton relaxation parameter lowered to %g <<<"
            print_text(s % par['relaxation_parameter'], 'red', 1)
            print_text(stars, 'red', 1)

      # solve velocity :
      else:
        momentum.solve(annotate=annotate)
    
      # solve mass equations, lowering time step on failure :
      if adaptive:
        solved_h = False
        while not solved_h:
          if dt < DOLFIN_EPS:
            status_h = [False,False]
            break
          H        = self.H.copy(True)
          status_h = mass.solve(annotate=annotate)
          solved_h = status_h[1]
          if t <= 100:
            solved_h = True
          if not solved_h:
            dt /= 2.0
            print_text(stars, 'red', 1)
            s = ">>> WARNING: time step lowered to %g <<<"
            print_text(s % dt, 'red', 1)
            self.init_time_step(dt, cls=self)
            self.init_H_H0(H, cls=self)
            print_text(stars, 'red', 1)

      # solve mass :
      else:
        mass.solve(annotate=annotate)
      
      ## use adaptive solver if desired :
      #if adaptive and (not mom_s[1] or not mas_s[1]):
      #  s = "::: reducing time step for convergence :::"
      #  print_text(s, self.color())
      #  solved, dt, t = self.adaptive_update(momentum, energy, mass,
      #                                       t_start, t_end, t,
      #                                       annotate=annotate)
      #  time_step = dt
      #  self.init_time_step(dt)
      
      # solve energy :
      energy.solve(annotate=annotate)

      # update pressure-melting point :
      energy.calc_T_melt(annotate=annotate)

      if callback != None:
        s    = '::: calling callback function :::'
        print_text(s, cls=self.this)
        callback()
       
      # increment time step :
      s = '>>> Time: %g yr, CPU time for last dt: %.3f s <<<'
      print_text(s % (t+dt, time()-tic), 'red', 1)

      t += dt
      self.step_time.append(time() - tic)
      
      # for the subsequent iteration, reset the parameters to normal :
      if adaptive:
        if par['relaxation_parameter'] != alpha:
          print_text("::: resetting alpha to normal :::", cls=self.this)
          par['relaxation_parameter'] = alpha
        if dt != time_step:
          print_text("::: resetting dt to normal :::", cls=self.this)
          self.init_time_step(time_step, cls=self.this)
          dt = time_step
      

    # calculate total time to compute
    s = time() - t0
    m = s / 60.0
    h = m / 60.0
    s = s % 60
    m = m % 60
    text = "total time to perform transient run: %02d:%02d:%02d" % (h,m,s)
    print_text(text, 'red', 1)



