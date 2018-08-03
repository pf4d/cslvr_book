from cslvr   import *

# set the relavent directories :
base_dir = 'dump/jakob_small/inversion_Wc_0.01/10/'
in_dir   = base_dir
out_dir  = base_dir + 'stress_balance/'
var_dir  = 'dump/vars_jakobshavn_small/'

# create HDF5 files for saving and loading data :
fdata   = HDF5File(mpi_comm_world(), var_dir + 'state.h5',       'r')
fin     = HDF5File(mpi_comm_world(), in_dir  + 'inverted_10.h5', 'r')
fout    = HDF5File(mpi_comm_world(), in_dir  + 'stress.h5',      'w')

# create 3D model for stokes solves :
model = D3Model(fdata, out_dir)

# init subdomains :
model.set_subdomains(fdata)

# initialize the 3D model vars :
model.init_S(fdata)
model.init_B(fdata)
model.init_mask(fdata)
model.init_adot(fdata)
model.init_beta(fin)
model.init_U(fin)
model.init_p(fin)
model.init_T(fin)
model.init_W(fin)
model.init_theta(fin)
model.form_energy_dependent_rate_factor()
model.calc_A()

mom = MomentumDukowiczStokes(model)
F   = StressBalance(model, momentum=mom)
F.solve()


model.save_hdf5(model.N_ii, f=fout)
model.save_hdf5(model.N_ij, f=fout)
model.save_hdf5(model.N_iz, f=fout)
model.save_hdf5(model.N_ji, f=fout)
model.save_hdf5(model.N_jj, f=fout)
model.save_hdf5(model.N_jz, f=fout)
model.save_hdf5(model.N_zi, f=fout)
model.save_hdf5(model.N_zj, f=fout)
model.save_hdf5(model.N_zz, f=fout)

model.save_hdf5(model.M_ii, f=fout)
model.save_hdf5(model.M_ij, f=fout)
model.save_hdf5(model.M_iz, f=fout)
model.save_hdf5(model.M_ji, f=fout)
model.save_hdf5(model.M_jj, f=fout)
model.save_hdf5(model.M_jz, f=fout)
model.save_hdf5(model.M_zi, f=fout)
model.save_hdf5(model.M_zj, f=fout)
model.save_hdf5(model.M_zz, f=fout)

fout.close()



