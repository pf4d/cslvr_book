from cslvr import *

alpha = 0.5 * pi / 180 
L     = 8000

p1    = Point(0.0, 0.0, 0.0)
p2    = Point(L,   L,   1)
mesh  = BoxMesh(p1, p2, 15, 15, 5)


out_dir = './ISMIP_HOM_A_results/RS/'
plt_dir = '../../../../images/momentum/ISMIP_HOM_A/small/RS/'

model = D3Model(mesh, out_dir = out_dir, use_periodic = True)

surface = Expression('- x[0] * tan(alpha)', alpha=alpha,
                     element=model.Q.ufl_element())
bed     = Expression(  '- x[0] * tan(alpha) - 1000.0 + 500.0 * ' \
                     + ' sin(2*pi*x[0]/L) * sin(2*pi*x[1]/L)',
                     alpha=alpha, L=L, element=model.Q.ufl_element())

model.calculate_boundaries()
model.deform_mesh_to_geometry(surface, bed)

model.init_mask(1.0)  # all grounded
model.init_beta(1000)
model.init_A(1e-16)
model.init_E(1.0)

#mom = MomentumDukowiczBP(model)
mom = MomentumDukowiczStokesReduced(model)
#mom = MomentumDukowiczStokes(model)
mom.solve()

divU = project(div(model.U3))

model.save_xdmf(model.p,  'p')
model.save_xdmf(model.U3, 'U')
model.save_xdmf(divU,     'divU')

#===============================================================================
# plotting :

# create the bed and surface meshes :
model.form_bed_mesh()
model.form_srf_mesh()

# create 2D models :
bedmodel = D2Model(model.bedmesh, out_dir)
srfmodel = D2Model(model.srfmesh, out_dir)

divU_b   = Function(bedmodel.Q)
bedmodel.assign_submesh_variable(divU_b, divU)

srfmodel.assign_submesh_variable(srfmodel.U3, model.U3)
srfmodel.init_U_mag(srfmodel.U3)
bedmodel.assign_submesh_variable(bedmodel.p,  model.p)

U_min  = srfmodel.U_mag.vector().min()
U_max  = srfmodel.U_mag.vector().max()
U_lvls = array([84, 86, 88, 90, 92, 94, 96, 98, 100])

p_min  = bedmodel.p.vector().min()
p_max  = bedmodel.p.vector().max()
p_lvls = array([4e6, 5e6, 6e6, 7e6, 8e6, 9e6, 1e7, 1.1e7, 1.2e7, p_max])

d_min  = divU_b.vector().min()
d_max  = divU_b.vector().max()
d_lvls = array([d_min, -5e-3, -2.5e-3, -1e-3, 
                1e-3, 2.5e-3, 5e-3, d_max])

plot_variable(u = srfmodel.U3, name = 'U_mag', direc = plt_dir,
              levels              = U_lvls,
              cmap                = 'viridis',
              tp                  = True,
              show                = False,
              extend              = 'both',#'neither',
              cb_format           = '%g')

plot_variable(u = bedmodel.p, name = 'p', direc = plt_dir,
              levels              = p_lvls,
              cmap                = 'viridis',
              tp                  = True,
              show                = False,
              extend              = 'min',
              cb_format           = '%.1e')

plot_variable(u = divU_b, name = 'divU', direc = plt_dir,
              cmap                = 'RdGy',
              levels              = d_lvls,
              tp                  = True,
              show                = False,
              extend              = 'neither',
              cb_format           = '%.1e')



