from cslvr import *

thklim = 1e-2
mesh_H = 5

# collect the raw data :
searise = DataFactory.get_searise(thklim)
bamber  = DataFactory.get_bamber(thklim)
rignot   = DataFactory.get_rignot()

# set the output directory :
out_dir = 'dump/vars/'

# load a mesh :
mesh  = Mesh('dump/meshes/greenland_2D_%iH_mesh.xml.gz' % mesh_H)

# create data objects to use with varglas :
dsr   = DataInput(searise, mesh=mesh)
dbm   = DataInput(bamber,  mesh=mesh)
drg   = DataInput(rignot,  mesh=mesh)

# Rignot dataset is on a different projection :
drg.change_projection(dbm)

B     = dbm.get_expression("B",     near=False)
S     = dbm.get_expression("S",     near=False)
M     = dbm.get_expression('mask',  near=True)
adot  = dsr.get_expression("adot",  near=False)
u_ob  = drg.get_expression("vx",    near=False)
v_ob  = drg.get_expression("vy",    near=False)
U_msk = drg.get_expression('mask',  near=True)

model = D2Model(mesh, out_dir = 'results/')
model.calculate_boundaries(mask=M, U_mask=U_msk, adot=adot) 

# calculate_boundaries() initializes model.mask, model.U_mask, and model.adot

model.init_S(S)
model.init_B(B)
model.init_U_ob(u_ob, v_ob)

lst = [model.S,
       model.B,
       model.adot,
       model.mask,
       model.u_ob,
       model.v_ob,
       model.U_mask]

f = HDF5File(mpi_comm_world(), out_dir + 'state.h5', 'w')

model.save_list_to_hdf5(lst, f)
model.save_subdomain_data(f)
model.save_mesh(f)

f.close()

model.save_xdmf(model.cf, 'cf')


