from cslvr             import *
from pylab             import *
from scipy.interpolate import RectBivariateSpline


#===============================================================================
# data preparation :
out_dir   = 'dump/meshes/'
mesh_name = 'jakobshavn_3D_small_block'

# get the data :
bamber   = DataFactory.get_bamber()
rignot   = DataFactory.get_rignot()

# process the data into the format used by CSLVR :
dbm      = DataInput(bamber,  gen_space=False)
drg      = DataInput(rignot,  gen_space=False)

# the `Rignot' dataset is defined with a different projection :
drg.change_projection(dbm)

# get surface velocity magnitude :
U_ob = sqrt(drg.data['vx']**2 + drg.data['vy']**2 + 1e-16)
drg.data['U_ob'] = U_ob


# form field from which to refine :
#===============================================================================

# this method will create a new data field with key 'ref' that we will 
# use to define the cell diameter, with maximum diameter umax, minimum
# diameter umin, and varied inversely proportional to the first argument;
# we want the cell size to be small where the ice is moving fast, and large
# where the ice is moving slowly :
drg.rescale_field('U_ob', 'ref', umin=500.0, umax=300000.0, inverse=True)

# eliminate just the edge of the mask so that we can properly interpolate
# the geometry to the terminus :
L = dbm.data['lat_mask']
dbm.data['mask'][L > 0.0] = 0


# generate the contour :
#===============================================================================

# the meshgenerator will create a mesh from the `Bamber' data :
m = MeshGenerator(dbm, mesh_name, out_dir)

# generate a contour around the land-ice mask `mask', skipping zero points :
m.create_contour('mask', zero_cntr=0.0001, skip_pts=0)

# coordinates of the box in x,y projection coordinates that we'll
# slice out of the continent :
x1 = -500000; y1 = -2190000
x2 = -150000; y2 = -2320000

# the x,y contour of the box `slice' :
new_cont = array([[x1, y1],
                  [x2, y1],
                  [x2, y2],
                  [x1, y2],
                  [x1, y1]])

# get the intersection of the contour of the entire continent and the box :
m.intersection(new_cont)

# make sure that no paths intersect of the new contour :
m.eliminate_intersections(dist=20)

# transform the contour to the projection used by the `Rignot' dataset :
m.transform_contour(rignot)

# make sure no two points lie directly on top of each other :
m.check_dist()

# save the resulting, perfect contour :
m.write_gmsh_contour(boundary_extend=False)

# you can plot the contour to make sure its perfect too :
#m.plot_contour()

# we have a 3D mesh here, so we extrude it any uniform distance vertically,
# with a number of layers :
m.extrude(h=100000, n_layers=10)

# close that .geo contour file so we don't corrupt any data :
m.close_file()


# mesh refinement :
#===============================================================================

# the MeshRefiner class uses a data array to set the cellsize to something
# we want, here the observed velocity speed of the `Rignot' dataset :
ref_bm = MeshRefiner(drg, 'ref', gmsh_file_name = out_dir + mesh_name)

# we need the background field from which we refine :
a,aid = ref_bm.add_static_attractor()
ref_bm.set_background_field(aid)

# this refines off of the background field we just set :
ref_bm.finish(gui = False, out_file_name = out_dir + mesh_name)

# and convert to the xml.gz file utilized by FEniCS :
ref_bm.convert_msh_to_xml()



