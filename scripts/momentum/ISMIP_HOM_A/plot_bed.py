from cslvr import *

# create a mesh and basal geometry identical to our problem, but now 
# in two dimension, cause we only want to plot the basal geometry pattern :
p1      = Point(0.0, 0.0)
p2      = Point(1.0, 1.0)
mesh    = RectangleMesh(p1, p2, 15, 15)

# where we're going to save that sweet plot :
plt_dir = '../../../images/momentum/ISMIP_HOM_A/'

# only need a 2D-model for this task, nevermind periodic boundaries too :
model   = D2Model(mesh, out_dir = plt_dir)
bed     = Expression('500.0 * sin(2*pi*x[0]/L) * sin(2*pi*x[1]/L)',
                     L=1, element=model.Q.ufl_element())

# initialize the one thing we care about (we could do other stuff too) :
model.init_B(bed)

# figure out the levels and plot them
B_min  = model.B.vector().min()  # we know this, but let's be sure
B_max  = model.B.vector().max()  # we know this, but let's be sure
B_lvls = array([B_min, -400, -300, -200, -100, -25,
                25, 100,  200,  300,  400,  B_max])

# this time, let's plot the topography like a topographic map :
plot_variable(u = model.B, name = 'B', direc = plt_dir,
              figsize = (5,5), levels = B_lvls, tp = True,
              show = False, cb = False, contour_type = 'lines',
              hide_ax_tick_labels = True)



