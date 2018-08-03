from cslvr   import *
from fenics  import *
from pylab   import *

out_dir  = 'dump/vars_jakobshavn_small/'

# collect the raw data :
searise  = DataFactory.get_searise()
bamber   = DataFactory.get_bamber(1.0)
rignot   = DataFactory.get_rignot()

# define the mesh :
mesh = Mesh('dump/meshes/jakobshavn_3D_small_block.xml.gz')

# create data objects to use with varglas :
dsr     = DataInput(searise,  mesh=mesh)
dbm     = DataInput(bamber,   mesh=mesh)
drg     = DataInput(rignot,   mesh=mesh)

# change the projection of all data to be the same as the mesh :
dbm.change_projection(drg)
dsr.change_projection(drg)

# get the expressions used by varglas :
S     = dbm.get_expression('S',        near=False)   # surface height
B     = dbm.get_expression('B',        near=False)   # bed height
M     = dbm.get_expression('mask',     near=True)    # shelf mask
L     = dbm.get_expression('lat_mask', near=True)    # terminus mask
adot  = dsr.get_expression('adot',     near=False)   # acc/abl function
q_geo = dsr.get_expression('q_geo',    near=False)   # geothermal heat
T_s   = dsr.get_expression('T',        near=False)   # surface temperature
u_ob  = drg.get_expression('vx',       near=False)   # x-comp. surface velocity
v_ob  = drg.get_expression('vy',       near=False)   # y-comp. surface velocity
U_msk = drg.get_expression('mask',     near=True)    # where U obs. are pres.

# create the 3D model :
model = D3Model(mesh=mesh, out_dir=out_dir)

# deform the mesh to match the geometry provided by the `Bamber' dataset :
model.deform_mesh_to_geometry(S, B)

# calculate the boundaries appropriately.  The lat_mask marks the exterior
# cliff edges of the mesh, U_mask marks areas without U observations, 
# adot creates the age equations surface boundary condition, and 
# mark_divide is set to True to in inform the model to mark the 
# boundaries interior to the domain.
model.calculate_boundaries(mask=M, lat_mask=L, U_mask=U_msk, adot=adot, 
                           mark_divide=True)

# initialize the observations :
model.init_T_surface(T_s)
model.init_q_geo(q_geo)
model.init_U_ob(u_ob, v_ob)

# these area all data fields we'll need to perform the data-assimilation
# procedure :
lst = [model.S,
       model.B,
       model.mask,
       model.q_geo,
       model.T_surface,
       model.adot,
       model.u_ob,
       model.v_ob,
       model.U_mask,
       model.lat_mask]

# the container for the data :
f = HDF5File(mpi_comm_world(), out_dir + 'state.h5', 'w')

# save all the data :
model.save_list_to_hdf5(lst, f)

# save the facet and cell markers generated by calculate_boundaries() :
model.save_subdomain_data(f)

# save the geometry-deformed mesh :
model.save_mesh(f)

f.close()


