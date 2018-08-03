from cslvr import *

# constants used :
a     = 0.5 * pi / 180 
L     = 8000

# create the mesh :
p1    = Point(0.0, 0.0, 0.0)
p2    = Point(L,   L,   1)
mesh  = BoxMesh(p1, p2, 15, 15, 5)

out_dir = './results/'

# stress balance requires the solution to a 3D model, in this case periodic :
model = D3Model(mesh, out_dir = out_dir, use_periodic = True)

# the topography expressions :
surface = Expression('- x[0] * tan(a)', a=a,
                     element=model.Q.ufl_element())
bed     = Expression(  '- x[0] * tan(a) - 1000.0 + 500.0 * ' \
                     + ' sin(2*pi*x[0]/L) * sin(2*pi*x[1]/L)',
                     a=a, L=L, element=model.Q.ufl_element())

# mark the boundaries used for the integration :
model.calculate_boundaries()

# deform the mesh to the topography :
model.deform_mesh_to_geometry(surface, bed)

model.init_beta(1000)   # traction
model.init_A(1e-16)     # isothermal rate-factor

# just have to change the momentum physics to balance the other stress models :
mom = MomentumDukowiczBP(model)
#mom = MomentumDukowiczStokesReduced(model)
#mom = MomentumDukowiczStokes(model)
mom.solve()

# solve the stress balance for the given Momentum instance :
F   = StressBalance(model, momentum=mom)
F.solve()

# save the resulting fields for plotting :
lst = [model.U3,
       model.p,
       model.N_ii,
       model.N_ij,
       model.N_iz,
       model.N_ji,
       model.N_jj,
       model.N_jz,
       model.N_zi,
       model.N_zj,
       model.N_zz,
       model.M_ii,
       model.M_ij,
       model.M_iz,
       model.M_ji,
       model.M_jj,
       model.M_jz,
       model.M_zi,
       model.M_zj,
       model.M_zz]

f = HDF5File(mpi_comm_world(), out_dir + 'FS.h5', 'w')

model.save_list_to_hdf5(lst, f)
model.save_subdomain_data(f)
model.save_mesh(f)

f.close()



