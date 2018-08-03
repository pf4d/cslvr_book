from fenics import *

mesh = Mesh("meshes/dolfin_fine.xml.gz")
sub_domains = MeshFunction("size_t", mesh, 
                           "meshes/dolfin_fine_subdomains.xml.gz")

# Taylor-Hood element
V     = VectorFunctionSpace(mesh, 'CG', 2)
Q     = FunctionSpace(mesh, 'CG', 1)
W     = V*Q

# variational problem :
u, p  = TrialFunctions(W)
v, q  = TestFunctions(W)

# no penetration boundary condition for velocity 
# y = 0, y = 1 and around the dolphin :
u_n   = Constant(0.0)

# inflow boundary condition for velocity at x = 1 :
u_0   = Expression(("-sin(x[1]*pi)", "0.0"))

# relavent measures :
ds    = Measure("ds")[sub_domains]
dG_0  = ds(0)
dG_r  = ds(1)

# constants :
gamma = Constant(1e2)
h     = CellSize(mesh)
n     = FacetNormal(mesh)
I     = Identity(2)
eta   = Constant(1.0)
f     = Constant((0.0,0.0))
beta  = Constant(10.0)

def epsilon(u): return 0.5*(grad(u) + grad(u).T)
def sigma(u,p): return 2*eta*epsilon(u) - p*I

t   = dot(sigma(u,p), n)
s   = dot(sigma(v,q), n)

B_o = + inner(sigma(u,p),grad(v))*dx - div(u)*q*dx

B_g = - dot(n,t) * dot(v,n) * dG_0 \
      - dot(u,n) * dot(s,n) * dG_0 \
      + gamma/h * dot(u,n) * dot(v,n) * dG_0 \
      + beta * dot(u, v) * dG_0 \
      - inner(dot(sigma(u,p), n), v) * dG_r \
      - inner(dot(sigma(v,q), n), u) * dG_r \
      + gamma/h * inner(v,u) * dG_r

F   = + dot(f,v) * dx \
      + gamma/h * u_n * dot(v,n) * dG_0 \
      - inner(dot(sigma(v,q), n), u_0) * dG_r \
      + gamma/h * inner(v,u_0) * dG_r

# solve variational problem
wh = Function(W)

solve(B_o + B_g == F, wh)
uh, ph = wh.split(True)

print "Norm of velocity coefficient vector: %.15g" % uh.vector().norm("l2")
print "Norm of pressure coefficient vector: %.15g" % ph.vector().norm("l2")

# get individual components with deep copy :
u0,u1 = uh.split(True)

from pylab                   import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib              import colors, ticker

# calculate array componets :
v0  = u0.compute_vertex_values(mesh)
v1  = u1.compute_vertex_values(mesh)
v   = sqrt(v0**2 + v1**2 + 1e-16)
v0  = v0 / v
v1  = v1 / v
x   = mesh.coordinates()[:,0]
y   = mesh.coordinates()[:,1]
t   = mesh.cells()

# generate velocity figure :
fig = figure(figsize=(8,7))
ax  = fig.add_subplot(111)

v[v > 2.0] = 2.0
cm = get_cmap('viridis')
ls = array([0.0,0.1,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.00001])
nm = colors.BoundaryNorm(ls, cm.N)
c  = ax.tricontourf(x, y, t, v, cmap=cm, norm=nm, levels=ls)
tp = ax.triplot(x, y, t, '-', color='k', lw=0.2, alpha=0.3)
q  = ax.quiver(x, y, v0, v1, pivot='middle',
                             color='k',
                             scale=60,
                             width=0.0015,
                             headwidth=4.0, 
                             headlength=4.0, 
                             headaxislength=4.0)
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.axis('equal')
ax.set_xlim([x.min(), x.max()])
ax.set_ylim([y.min(), y.max()])
ax.set_xticklabels([])
ax.set_yticklabels([])
  
divider = make_axes_locatable(gca())
cax  = divider.append_axes('right', "5%", pad="3%")
cbar = fig.colorbar(c, cax=cax, ticks=ls, format='%.1f') 
tight_layout()
savefig('../../images/fenics_intro/2Dstokes_nitsche_u.pdf')

# generate pressure figure :
v  = ph.compute_vertex_values(mesh)

fig = figure(figsize=(8,7))
ax  = fig.add_subplot(111)

v[v >  120] = 120
v[v < -20] = -20

ls = array([v.min(),-10, 0,10,20,30,40,50,60,70,80,90,100,110,120])
nm = colors.BoundaryNorm(ls, cm.N)
c  = ax.tricontourf(x, y, t, v, 10, cmap=cm, norm=nm, levels=ls)
tp = ax.triplot(x, y, t, '-', color='k', lw=0.2, alpha=0.5)

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.axis('equal')
ax.set_xlim([x.min(), x.max()])
ax.set_ylim([y.min(), y.max()])
ax.set_xticklabels([])
ax.set_yticklabels([])
  
divider = make_axes_locatable(gca())
cax  = divider.append_axes('right', "5%", pad="3%")
cbar = colorbar(c, cax=cax) 
tight_layout()
savefig('../../images/fenics_intro/2Dstokes_nitsche_p.pdf')



