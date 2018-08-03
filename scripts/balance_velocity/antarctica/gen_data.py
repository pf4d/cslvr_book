from cslvr import *

thklim = 1.0
mesh_H = 10

# collect the raw data :
bedmap2  = DataFactory.get_bedmap2(thklim=thklim)
bedmap1  = DataFactory.get_bedmap1(thklim=thklim)
measures = DataFactory.get_ant_measures(res=900)

# set the output directory :
out_dir = 'dump/vars/'

# load a mesh :
mesh  = Mesh('dump/meshes/antarctica_2D_%iH_mesh.xml.gz' % mesh_H)

# create data objects to use with varglas :
db1   = DataInput(bedmap1,  mesh=mesh)
db2   = DataInput(bedmap2,  mesh=mesh)
dbm   = DataInput(measures, mesh=mesh)

S     = db2.get_expression("S",        near=False)
B     = db2.get_expression("B",        near=False)
M     = db2.get_expression("mask",     near=True)
adot  = db1.get_expression("acca",     near=False)
u_ob  = dbm.get_expression("vx",       near=False)
v_ob  = dbm.get_expression("vy",       near=False)
U_msk = dbm.get_expression("mask",     near=True)

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



