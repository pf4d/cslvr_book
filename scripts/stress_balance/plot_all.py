from cslvr import *

fs_out_dir = '../../images/stress_balance/FS/'
rs_out_dir = '../../images/stress_balance/RS/'
bp_out_dir = '../../images/stress_balance/BP/'
in_dir     = './results/'

ffs  = HDF5File(mpi_comm_world(), in_dir + 'FS.h5', 'r')
frs  = HDF5File(mpi_comm_world(), in_dir + 'RS.h5', 'r')
fbp  = HDF5File(mpi_comm_world(), in_dir + 'BP.h5', 'r')

fs_model = D3Model(ffs, out_dir = fs_out_dir)
rs_model = D3Model(frs, out_dir = rs_out_dir)
bp_model = D3Model(fbp, out_dir = bp_out_dir)

fs_model.set_subdomains(ffs)
rs_model.set_subdomains(frs)
bp_model.set_subdomains(fbp)

fs_model.init_U(ffs)
fs_model.init_p(ffs)
fs_model.init_N_ii(ffs)
fs_model.init_N_ij(ffs)
fs_model.init_N_iz(ffs)
fs_model.init_N_ji(ffs)
fs_model.init_N_jj(ffs)
fs_model.init_N_jz(ffs)
fs_model.init_N_zi(ffs)
fs_model.init_N_zj(ffs)
fs_model.init_N_zz(ffs)
fs_model.init_M_ii(ffs)
fs_model.init_M_ij(ffs)
fs_model.init_M_iz(ffs)
fs_model.init_M_ji(ffs)
fs_model.init_M_jj(ffs)
fs_model.init_M_jz(ffs)
fs_model.init_M_zi(ffs)
fs_model.init_M_zj(ffs)
fs_model.init_M_zz(ffs)

rs_model.init_U(frs)
rs_model.init_p(frs)
rs_model.init_N_ii(frs)
rs_model.init_N_ij(frs)
rs_model.init_N_iz(frs)
rs_model.init_N_ji(frs)
rs_model.init_N_jj(frs)
rs_model.init_N_jz(frs)
rs_model.init_N_zi(frs)
rs_model.init_N_zj(frs)
rs_model.init_N_zz(frs)
rs_model.init_M_ii(frs)
rs_model.init_M_ij(frs)
rs_model.init_M_iz(frs)
rs_model.init_M_ji(frs)
rs_model.init_M_jj(frs)
rs_model.init_M_jz(frs)
rs_model.init_M_zi(frs)
rs_model.init_M_zj(frs)
rs_model.init_M_zz(frs)

bp_model.init_U(fbp)
bp_model.init_p(fbp)
bp_model.init_N_ii(fbp)
bp_model.init_N_ij(fbp)
bp_model.init_N_iz(fbp)
bp_model.init_N_ji(fbp)
bp_model.init_N_jj(fbp)
bp_model.init_N_jz(fbp)
bp_model.init_N_zi(fbp)
bp_model.init_N_zj(fbp)
bp_model.init_N_zz(fbp)
bp_model.init_M_ii(fbp)
bp_model.init_M_ij(fbp)
bp_model.init_M_iz(fbp)
bp_model.init_M_ji(fbp)
bp_model.init_M_jj(fbp)
bp_model.init_M_jz(fbp)
bp_model.init_M_zi(fbp)
bp_model.init_M_zj(fbp)
bp_model.init_M_zz(fbp)

#===============================================================================
# plotting :

# create the bed and surface meshes :
fs_model.form_bed_mesh()
fs_model.form_srf_mesh()
rs_model.form_bed_mesh()
rs_model.form_srf_mesh()
bp_model.form_bed_mesh()
bp_model.form_srf_mesh()

# create 2D models :
fs_bedmodel = D2Model(fs_model.bedmesh, fs_out_dir)
fs_srfmodel = D2Model(fs_model.srfmesh, fs_out_dir)
rs_bedmodel = D2Model(rs_model.bedmesh, rs_out_dir)
rs_srfmodel = D2Model(rs_model.srfmesh, rs_out_dir)
bp_bedmodel = D2Model(bp_model.bedmesh, bp_out_dir)
bp_srfmodel = D2Model(bp_model.srfmesh, bp_out_dir)

fs_srfmodel.assign_submesh_variable(fs_srfmodel.U_mag,  fs_model.U_mag)
fs_srfmodel.assign_submesh_variable(fs_srfmodel.U3,     fs_model.U3)
fs_bedmodel.assign_submesh_variable(fs_bedmodel.p,      fs_model.p)
fs_srfmodel.assign_submesh_variable(fs_srfmodel.N_ii,   fs_model.N_ii)
fs_srfmodel.assign_submesh_variable(fs_srfmodel.N_ij,   fs_model.N_ij)
fs_srfmodel.assign_submesh_variable(fs_srfmodel.N_iz,   fs_model.N_iz)
fs_srfmodel.assign_submesh_variable(fs_srfmodel.N_ji,   fs_model.N_ji)
fs_srfmodel.assign_submesh_variable(fs_srfmodel.N_jj,   fs_model.N_jj)
fs_srfmodel.assign_submesh_variable(fs_srfmodel.N_jz,   fs_model.N_jz)
fs_srfmodel.assign_submesh_variable(fs_srfmodel.N_zi,   fs_model.N_zi)
fs_srfmodel.assign_submesh_variable(fs_srfmodel.N_zj,   fs_model.N_zj)
fs_srfmodel.assign_submesh_variable(fs_srfmodel.N_zz,   fs_model.N_zz)
fs_srfmodel.assign_submesh_variable(fs_srfmodel.M_ii,   fs_model.M_ii)
fs_srfmodel.assign_submesh_variable(fs_srfmodel.M_ij,   fs_model.M_ij)
fs_srfmodel.assign_submesh_variable(fs_srfmodel.M_iz,   fs_model.M_iz)
fs_srfmodel.assign_submesh_variable(fs_srfmodel.M_ji,   fs_model.M_ji)
fs_srfmodel.assign_submesh_variable(fs_srfmodel.M_jj,   fs_model.M_jj)
fs_srfmodel.assign_submesh_variable(fs_srfmodel.M_jz,   fs_model.M_jz)
fs_srfmodel.assign_submesh_variable(fs_srfmodel.M_zi,   fs_model.M_zi)
fs_srfmodel.assign_submesh_variable(fs_srfmodel.M_zj,   fs_model.M_zj)
fs_srfmodel.assign_submesh_variable(fs_srfmodel.M_zz,   fs_model.M_zz)

rs_srfmodel.assign_submesh_variable(rs_srfmodel.U_mag,  rs_model.U_mag)
rs_srfmodel.assign_submesh_variable(rs_srfmodel.U3,     rs_model.U3)
rs_bedmodel.assign_submesh_variable(rs_bedmodel.p,      rs_model.p)
rs_srfmodel.assign_submesh_variable(rs_srfmodel.N_ii,   rs_model.N_ii)
rs_srfmodel.assign_submesh_variable(rs_srfmodel.N_ij,   rs_model.N_ij)
rs_srfmodel.assign_submesh_variable(rs_srfmodel.N_iz,   rs_model.N_iz)
rs_srfmodel.assign_submesh_variable(rs_srfmodel.N_ji,   rs_model.N_ji)
rs_srfmodel.assign_submesh_variable(rs_srfmodel.N_jj,   rs_model.N_jj)
rs_srfmodel.assign_submesh_variable(rs_srfmodel.N_jz,   rs_model.N_jz)
rs_srfmodel.assign_submesh_variable(rs_srfmodel.N_zi,   rs_model.N_zi)
rs_srfmodel.assign_submesh_variable(rs_srfmodel.N_zj,   rs_model.N_zj)
rs_srfmodel.assign_submesh_variable(rs_srfmodel.N_zz,   rs_model.N_zz)
rs_srfmodel.assign_submesh_variable(rs_srfmodel.M_ii,   rs_model.M_ii)
rs_srfmodel.assign_submesh_variable(rs_srfmodel.M_ij,   rs_model.M_ij)
rs_srfmodel.assign_submesh_variable(rs_srfmodel.M_iz,   rs_model.M_iz)
rs_srfmodel.assign_submesh_variable(rs_srfmodel.M_ji,   rs_model.M_ji)
rs_srfmodel.assign_submesh_variable(rs_srfmodel.M_jj,   rs_model.M_jj)
rs_srfmodel.assign_submesh_variable(rs_srfmodel.M_jz,   rs_model.M_jz)
rs_srfmodel.assign_submesh_variable(rs_srfmodel.M_zi,   rs_model.M_zi)
rs_srfmodel.assign_submesh_variable(rs_srfmodel.M_zj,   rs_model.M_zj)
rs_srfmodel.assign_submesh_variable(rs_srfmodel.M_zz,   rs_model.M_zz)

bp_srfmodel.assign_submesh_variable(bp_srfmodel.U_mag,  bp_model.U_mag)
bp_srfmodel.assign_submesh_variable(bp_srfmodel.U3,     bp_model.U3)
bp_bedmodel.assign_submesh_variable(bp_bedmodel.p,      bp_model.p)
bp_srfmodel.assign_submesh_variable(bp_srfmodel.N_ii,   bp_model.N_ii)
bp_srfmodel.assign_submesh_variable(bp_srfmodel.N_ij,   bp_model.N_ij)
bp_srfmodel.assign_submesh_variable(bp_srfmodel.N_iz,   bp_model.N_iz)
bp_srfmodel.assign_submesh_variable(bp_srfmodel.N_ji,   bp_model.N_ji)
bp_srfmodel.assign_submesh_variable(bp_srfmodel.N_jj,   bp_model.N_jj)
bp_srfmodel.assign_submesh_variable(bp_srfmodel.N_jz,   bp_model.N_jz)
bp_srfmodel.assign_submesh_variable(bp_srfmodel.N_zi,   bp_model.N_zi)
bp_srfmodel.assign_submesh_variable(bp_srfmodel.N_zj,   bp_model.N_zj)
bp_srfmodel.assign_submesh_variable(bp_srfmodel.N_zz,   bp_model.N_zz)
bp_srfmodel.assign_submesh_variable(bp_srfmodel.M_ii,   bp_model.M_ii)
bp_srfmodel.assign_submesh_variable(bp_srfmodel.M_ij,   bp_model.M_ij)
bp_srfmodel.assign_submesh_variable(bp_srfmodel.M_iz,   bp_model.M_iz)
bp_srfmodel.assign_submesh_variable(bp_srfmodel.M_ji,   bp_model.M_ji)
bp_srfmodel.assign_submesh_variable(bp_srfmodel.M_jj,   bp_model.M_jj)
bp_srfmodel.assign_submesh_variable(bp_srfmodel.M_jz,   bp_model.M_jz)
bp_srfmodel.assign_submesh_variable(bp_srfmodel.M_zi,   bp_model.M_zi)
bp_srfmodel.assign_submesh_variable(bp_srfmodel.M_zj,   bp_model.M_zj)
bp_srfmodel.assign_submesh_variable(bp_srfmodel.M_zz,   bp_model.M_zz)

fs_U_min    = fs_srfmodel.U_mag.vector().min()
fs_U_max    = fs_srfmodel.U_mag.vector().max()
fs_p_min    = fs_bedmodel.p.vector().min()
fs_p_max    = fs_bedmodel.p.vector().max()
fs_N_ii_min = fs_srfmodel.N_ii.vector().min()
fs_N_ii_max = fs_srfmodel.N_ii.vector().max()
fs_N_ij_min = fs_srfmodel.N_ij.vector().min()
fs_N_ij_max = fs_srfmodel.N_ij.vector().max()
fs_N_iz_min = fs_srfmodel.N_iz.vector().min()
fs_N_iz_max = fs_srfmodel.N_iz.vector().max()
fs_N_ji_min = fs_srfmodel.N_ji.vector().min()
fs_N_ji_max = fs_srfmodel.N_ji.vector().max()
fs_N_jj_min = fs_srfmodel.N_jj.vector().min()
fs_N_jj_max = fs_srfmodel.N_jj.vector().max()
fs_N_jz_min = fs_srfmodel.N_jz.vector().min()
fs_N_jz_max = fs_srfmodel.N_jz.vector().max()
fs_N_zi_min = fs_srfmodel.N_zi.vector().min()
fs_N_zi_max = fs_srfmodel.N_zi.vector().max()
fs_N_zj_min = fs_srfmodel.N_zj.vector().min()
fs_N_zj_max = fs_srfmodel.N_zj.vector().max()
fs_N_zz_min = fs_srfmodel.N_zz.vector().min()
fs_N_zz_max = fs_srfmodel.N_zz.vector().max()
fs_M_ii_min = fs_srfmodel.M_ii.vector().min()
fs_M_ii_max = fs_srfmodel.M_ii.vector().max()
fs_M_ij_min = fs_srfmodel.M_ij.vector().min()
fs_M_ij_max = fs_srfmodel.M_ij.vector().max()
fs_M_iz_min = fs_srfmodel.M_iz.vector().min()
fs_M_iz_max = fs_srfmodel.M_iz.vector().max()
fs_M_ji_min = fs_srfmodel.M_ji.vector().min()
fs_M_ji_max = fs_srfmodel.M_ji.vector().max()
fs_M_jj_min = fs_srfmodel.M_jj.vector().min()
fs_M_jj_max = fs_srfmodel.M_jj.vector().max()
fs_M_jz_min = fs_srfmodel.M_jz.vector().min()
fs_M_jz_max = fs_srfmodel.M_jz.vector().max()
fs_M_zi_min = fs_srfmodel.M_zi.vector().min()
fs_M_zi_max = fs_srfmodel.M_zi.vector().max()
fs_M_zj_min = fs_srfmodel.M_zj.vector().min()
fs_M_zj_max = fs_srfmodel.M_zj.vector().max()
fs_M_zz_min = fs_srfmodel.M_zz.vector().min()
fs_M_zz_max = fs_srfmodel.M_zz.vector().max()

rs_U_min    = rs_srfmodel.U_mag.vector().min()
rs_U_max    = rs_srfmodel.U_mag.vector().max()
rs_p_min    = rs_bedmodel.p.vector().min()
rs_p_max    = rs_bedmodel.p.vector().max()
rs_N_ii_min = rs_srfmodel.N_ii.vector().min()
rs_N_ii_max = rs_srfmodel.N_ii.vector().max()
rs_N_ij_min = rs_srfmodel.N_ij.vector().min()
rs_N_ij_max = rs_srfmodel.N_ij.vector().max()
rs_N_iz_min = rs_srfmodel.N_iz.vector().min()
rs_N_iz_max = rs_srfmodel.N_iz.vector().max()
rs_N_ji_min = rs_srfmodel.N_ji.vector().min()
rs_N_ji_max = rs_srfmodel.N_ji.vector().max()
rs_N_jj_min = rs_srfmodel.N_jj.vector().min()
rs_N_jj_max = rs_srfmodel.N_jj.vector().max()
rs_N_jz_min = rs_srfmodel.N_jz.vector().min()
rs_N_jz_max = rs_srfmodel.N_jz.vector().max()
rs_N_zi_min = rs_srfmodel.N_zi.vector().min()
rs_N_zi_max = rs_srfmodel.N_zi.vector().max()
rs_N_zj_min = rs_srfmodel.N_zj.vector().min()
rs_N_zj_max = rs_srfmodel.N_zj.vector().max()
rs_N_zz_min = rs_srfmodel.N_zz.vector().min()
rs_N_zz_max = rs_srfmodel.N_zz.vector().max()
rs_M_ii_min = rs_srfmodel.M_ii.vector().min()
rs_M_ii_max = rs_srfmodel.M_ii.vector().max()
rs_M_ij_min = rs_srfmodel.M_ij.vector().min()
rs_M_ij_max = rs_srfmodel.M_ij.vector().max()
rs_M_iz_min = rs_srfmodel.M_iz.vector().min()
rs_M_iz_max = rs_srfmodel.M_iz.vector().max()
rs_M_ji_min = rs_srfmodel.M_ji.vector().min()
rs_M_ji_max = rs_srfmodel.M_ji.vector().max()
rs_M_jj_min = rs_srfmodel.M_jj.vector().min()
rs_M_jj_max = rs_srfmodel.M_jj.vector().max()
rs_M_jz_min = rs_srfmodel.M_jz.vector().min()
rs_M_jz_max = rs_srfmodel.M_jz.vector().max()
rs_M_zi_min = rs_srfmodel.M_zi.vector().min()
rs_M_zi_max = rs_srfmodel.M_zi.vector().max()
rs_M_zj_min = rs_srfmodel.M_zj.vector().min()
rs_M_zj_max = rs_srfmodel.M_zj.vector().max()
rs_M_zz_min = rs_srfmodel.M_zz.vector().min()
rs_M_zz_max = rs_srfmodel.M_zz.vector().max()

bp_U_min    = bp_srfmodel.U_mag.vector().min()
bp_U_max    = bp_srfmodel.U_mag.vector().max()
bp_p_min    = bp_bedmodel.p.vector().min()
bp_p_max    = bp_bedmodel.p.vector().max()
bp_N_ii_min = bp_srfmodel.N_ii.vector().min()
bp_N_ii_max = bp_srfmodel.N_ii.vector().max()
bp_N_ij_min = bp_srfmodel.N_ij.vector().min()
bp_N_ij_max = bp_srfmodel.N_ij.vector().max()
bp_N_iz_min = bp_srfmodel.N_iz.vector().min()
bp_N_iz_max = bp_srfmodel.N_iz.vector().max()
bp_N_ji_min = bp_srfmodel.N_ji.vector().min()
bp_N_ji_max = bp_srfmodel.N_ji.vector().max()
bp_N_jj_min = bp_srfmodel.N_jj.vector().min()
bp_N_jj_max = bp_srfmodel.N_jj.vector().max()
bp_N_jz_min = bp_srfmodel.N_jz.vector().min()
bp_N_jz_max = bp_srfmodel.N_jz.vector().max()
bp_N_zi_min = bp_srfmodel.N_zi.vector().min()
bp_N_zi_max = bp_srfmodel.N_zi.vector().max()
bp_N_zj_min = bp_srfmodel.N_zj.vector().min()
bp_N_zj_max = bp_srfmodel.N_zj.vector().max()
bp_N_zz_min = bp_srfmodel.N_zz.vector().min()
bp_N_zz_max = bp_srfmodel.N_zz.vector().max()
bp_M_ii_min = bp_srfmodel.M_ii.vector().min()
bp_M_ii_max = bp_srfmodel.M_ii.vector().max()
bp_M_ij_min = bp_srfmodel.M_ij.vector().min()
bp_M_ij_max = bp_srfmodel.M_ij.vector().max()
bp_M_iz_min = bp_srfmodel.M_iz.vector().min()
bp_M_iz_max = bp_srfmodel.M_iz.vector().max()
bp_M_ji_min = bp_srfmodel.M_ji.vector().min()
bp_M_ji_max = bp_srfmodel.M_ji.vector().max()
bp_M_jj_min = bp_srfmodel.M_jj.vector().min()
bp_M_jj_max = bp_srfmodel.M_jj.vector().max()
bp_M_jz_min = bp_srfmodel.M_jz.vector().min()
bp_M_jz_max = bp_srfmodel.M_jz.vector().max()
bp_M_zi_min = bp_srfmodel.M_zi.vector().min()
bp_M_zi_max = bp_srfmodel.M_zi.vector().max()
bp_M_zj_min = bp_srfmodel.M_zj.vector().min()
bp_M_zj_max = bp_srfmodel.M_zj.vector().max()
bp_M_zz_min = bp_srfmodel.M_zz.vector().min()
bp_M_zz_max = bp_srfmodel.M_zz.vector().max()

bp_U_min    = bp_srfmodel.U_mag.vector().min()
bp_U_max    = bp_srfmodel.U_mag.vector().max()
bp_p_min    = bp_bedmodel.p.vector().min()
bp_p_max    = bp_bedmodel.p.vector().max()
bp_N_ii_min = bp_srfmodel.N_ii.vector().min()
bp_N_ii_max = bp_srfmodel.N_ii.vector().max()
bp_N_ij_min = bp_srfmodel.N_ij.vector().min()
bp_N_ij_max = bp_srfmodel.N_ij.vector().max()
bp_N_iz_min = bp_srfmodel.N_iz.vector().min()
bp_N_iz_max = bp_srfmodel.N_iz.vector().max()
bp_N_ji_min = bp_srfmodel.N_ji.vector().min()
bp_N_ji_max = bp_srfmodel.N_ji.vector().max()
bp_N_jj_min = bp_srfmodel.N_jj.vector().min()
bp_N_jj_max = bp_srfmodel.N_jj.vector().max()
bp_N_jz_min = bp_srfmodel.N_jz.vector().min()
bp_N_jz_max = bp_srfmodel.N_jz.vector().max()
bp_N_zi_min = bp_srfmodel.N_zi.vector().min()
bp_N_zi_max = bp_srfmodel.N_zi.vector().max()
bp_N_zj_min = bp_srfmodel.N_zj.vector().min()
bp_N_zj_max = bp_srfmodel.N_zj.vector().max()
bp_N_zz_min = bp_srfmodel.N_zz.vector().min()
bp_N_zz_max = bp_srfmodel.N_zz.vector().max()
bp_M_ii_min = bp_srfmodel.M_ii.vector().min()
bp_M_ii_max = bp_srfmodel.M_ii.vector().max()
bp_M_ij_min = bp_srfmodel.M_ij.vector().min()
bp_M_ij_max = bp_srfmodel.M_ij.vector().max()
bp_M_iz_min = bp_srfmodel.M_iz.vector().min()
bp_M_iz_max = bp_srfmodel.M_iz.vector().max()
bp_M_ji_min = bp_srfmodel.M_ji.vector().min()
bp_M_ji_max = bp_srfmodel.M_ji.vector().max()
bp_M_jj_min = bp_srfmodel.M_jj.vector().min()
bp_M_jj_max = bp_srfmodel.M_jj.vector().max()
bp_M_jz_min = bp_srfmodel.M_jz.vector().min()
bp_M_jz_max = bp_srfmodel.M_jz.vector().max()
bp_M_zi_min = bp_srfmodel.M_zi.vector().min()
bp_M_zi_max = bp_srfmodel.M_zi.vector().max()
bp_M_zj_min = bp_srfmodel.M_zj.vector().min()
bp_M_zj_max = bp_srfmodel.M_zj.vector().max()
bp_M_zz_min = bp_srfmodel.M_zz.vector().min()
bp_M_zz_max = bp_srfmodel.M_zz.vector().max()

U_min    = min(fs_U_min, rs_U_min, bp_U_min)
U_max    = max(fs_U_max, rs_U_max, bp_U_max)
p_min    = 4e6#min(fs_p_min, rs_p_min, bp_p_min)
p_max    = max(fs_p_max, rs_p_max, bp_p_max)
N_ii_min = min(fs_N_ii_min, rs_N_ii_min, bp_N_ii_min)
N_ii_max = max(fs_N_ii_max, rs_N_ii_max, bp_N_ii_max)
N_ij_min = min(fs_N_ij_min, rs_N_ij_min, bp_N_ij_min)
N_ij_max = max(fs_N_ij_max, rs_N_ij_max, bp_N_ij_max)
N_iz_min = min(fs_N_iz_min, rs_N_iz_min, bp_N_iz_min)
N_iz_max = max(fs_N_iz_max, rs_N_iz_max, bp_N_iz_max)
N_ji_min = min(fs_N_ji_min, rs_N_ji_min, bp_N_ji_min)
N_ji_max = max(fs_N_ji_max, rs_N_ji_max, bp_N_ji_max)
N_jj_min = min(fs_N_jj_min, rs_N_jj_min, bp_N_jj_min)
N_jj_max = max(fs_N_jj_max, rs_N_jj_max, bp_N_jj_max)
N_jz_min = min(fs_N_jz_min, rs_N_jz_min, bp_N_jz_min)
N_jz_max = max(fs_N_jz_max, rs_N_jz_max, bp_N_jz_max)
N_zi_min = min(fs_N_zi_min, rs_N_zi_min, bp_N_zi_min)
N_zi_max = max(fs_N_zi_max, rs_N_zi_max, bp_N_zi_max)
N_zj_min = min(fs_N_zj_min, rs_N_zj_min, bp_N_zj_min)
N_zj_max = max(fs_N_zj_max, rs_N_zj_max, bp_N_zj_max)
N_zz_min = min(fs_N_zz_min, rs_N_zz_min, bp_N_zz_min)
N_zz_max = max(fs_N_zz_max, rs_N_zz_max, bp_N_zz_max)
M_ii_min = min(fs_M_ii_min, rs_M_ii_min, bp_M_ii_min)
M_ii_max = max(fs_M_ii_max, rs_M_ii_max, bp_M_ii_max)
M_ij_min = min(fs_M_ij_min, rs_M_ij_min, bp_M_ij_min)
M_ij_max = max(fs_M_ij_max, rs_M_ij_max, bp_M_ij_max)
M_iz_min = min(fs_M_iz_min, rs_M_iz_min, bp_M_iz_min)
M_iz_max = max(fs_M_iz_max, rs_M_iz_max, bp_M_iz_max) * 1.5
M_ji_min = min(fs_M_ji_min, rs_M_ji_min, bp_M_ji_min)
M_ji_max = max(fs_M_ji_max, rs_M_ji_max, bp_M_ji_max)
M_jj_min = min(fs_M_jj_min, rs_M_jj_min, bp_M_jj_min)
M_jj_max = max(fs_M_jj_max, rs_M_jj_max, bp_M_jj_max)
M_jz_min = min(fs_M_jz_min, rs_M_jz_min, bp_M_jz_min)
M_jz_max = max(fs_M_jz_max, rs_M_jz_max, bp_M_jz_max)
M_zi_min = min(fs_M_zi_min, rs_M_zi_min, bp_M_zi_min)
M_zi_max = max(fs_M_zi_max, rs_M_zi_max, bp_M_zi_max)
M_zj_min = min(fs_M_zj_min, rs_M_zj_min, bp_M_zj_min)
M_zj_max = max(fs_M_zj_max, rs_M_zj_max, bp_M_zj_max)
M_zz_min = min(fs_M_zz_min, rs_M_zz_min, bp_M_zz_min)
M_zz_max = max(fs_M_zz_max, rs_M_zz_max, bp_M_zz_max)

U_a_max    = max(abs(U_min),    abs(U_max))
p_a_max    = max(abs(p_max),    abs(p_min))
N_ii_a_max = max(abs(N_ii_min), abs(N_ii_max))
N_ij_a_max = max(abs(N_ij_min), abs(N_ij_max))
N_iz_a_max = max(abs(N_iz_min), abs(N_iz_max))
N_ji_a_max = max(abs(N_ji_min), abs(N_ji_max))
N_jj_a_max = max(abs(N_jj_min), abs(N_jj_max))
N_jz_a_max = max(abs(N_jz_min), abs(N_jz_max))
N_zi_a_max = max(abs(N_zi_min), abs(N_zi_max))
N_zj_a_max = max(abs(N_zj_min), abs(N_zj_max))
N_zz_a_max = max(abs(N_zz_min), abs(N_zz_max))
M_ii_a_max = max(abs(M_ii_min), abs(M_ii_max)) / 1.5
M_ij_a_max = max(abs(M_ij_min), abs(M_ij_max)) / 1.5
M_iz_a_max = max(abs(M_iz_min), abs(M_iz_max)) / 2.0
M_ji_a_max = max(abs(M_ji_min), abs(M_ji_max)) / 3.0
M_jj_a_max = max(abs(M_jj_min), abs(M_jj_max))
M_jz_a_max = max(abs(M_jz_min), abs(M_jz_max)) / 2.0
M_zi_a_max = max(abs(M_zi_min), abs(M_zi_max))
M_zj_a_max = max(abs(M_zj_min), abs(M_zj_max)) / 2.0
M_zz_a_max = max(abs(M_zz_min), abs(M_zz_max))

N_ii_min = -N_ii_a_max
N_ii_max =  N_ii_a_max
N_ij_min = -N_ij_a_max
N_ij_max =  N_ij_a_max
#N_iz_min =  N_iz_a_max
#N_iz_max =  N_iz_a_max
N_ji_min = -N_ji_a_max
N_ji_max =  N_ji_a_max
N_jj_min = -N_jj_a_max
N_jj_max =  N_jj_a_max
N_jz_min = -N_jz_a_max
N_jz_max =  N_jz_a_max
#N_zi_min = -N_zi_a_max
#N_zi_max =  N_zi_a_max
N_zj_min = -N_zj_a_max
N_zj_max =  N_zj_a_max
N_zz_min = -N_zz_a_max
N_zz_max =  N_zz_a_max
M_ii_min = -M_ii_a_max
M_ii_max =  M_ii_a_max
M_ij_min = -M_ij_a_max
M_ij_max =  M_ij_a_max
#M_iz_min =  M_iz_a_max
#M_iz_max =  M_iz_a_max
M_ji_min = -M_ji_a_max
M_ji_max =  M_ji_a_max
M_jj_min = -M_jj_a_max
M_jj_max =  M_jj_a_max
M_jz_min = -M_jz_a_max
M_jz_max =  M_jz_a_max
M_zi_min = -M_zi_a_max
M_zi_max =  M_zi_a_max
M_zj_min = -M_zj_a_max
M_zj_max =  M_zj_a_max
M_zz_min = -M_zz_a_max
M_zz_max =  M_zz_a_max


#===============================================================================
# plot :

cmap     = 'viridis'
p_cmap   = 'gist_yarg'
n_cmap   = 'Reds_r'
sym_cmap = 'RdGy'

plot_variable(u                   = fs_srfmodel.U3,
              name                = 'U_mag',
              umin                = U_min,
              umax                = U_max,
              vec_scale           = 22,
              direc               = fs_out_dir,
              cmap                = cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              cb_format           = '%g')

plot_variable(u                   = fs_bedmodel.p,
              name                = 'p',
              umin                = p_min,
              umax                = p_max,
              direc               = fs_out_dir,
              cmap                = cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              cb_format           = '%.1e')

plot_variable(u                   = fs_srfmodel.N_ii,
              name                = 'N_ii',
              direc               = fs_out_dir,
              umin                = N_ii_min,
              umax                = N_ii_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

plot_variable(u                   = fs_srfmodel.N_ij,
              name                = 'N_ij',
              direc               = fs_out_dir,
              umin                = N_ij_min,
              umax                = N_ij_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

plot_variable(u                   = fs_srfmodel.N_iz,
              name                = 'N_iz',
              direc               = fs_out_dir,
              umin                = N_iz_min,
              umax                = N_iz_max,
              cmap                = p_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

plot_variable(u                   = fs_srfmodel.N_ji,
              name                = 'N_ji',
              direc               = fs_out_dir,
              umin                = N_ji_min,
              umax                = N_ji_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

plot_variable(u                   = fs_srfmodel.N_jj,
              name                = 'N_jj',
              direc               = fs_out_dir,
              umin                = N_jj_min,
              umax                = N_jj_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

plot_variable(u                   = fs_srfmodel.N_jz,
              name                = 'N_jz',
              direc               = fs_out_dir,
              umin                = N_jz_min,
              umax                = N_jz_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = fs_srfmodel.N_zi,
              name                = 'N_zi',
              direc               = fs_out_dir,
              umin                = N_zi_min,
              umax                = N_zi_max,
              cmap                = p_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = fs_srfmodel.N_zj,
              name                = 'N_zj',
              direc               = fs_out_dir,
              umin                = N_zj_min,
              umax                = N_zj_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = fs_srfmodel.N_zz,
              name                = 'N_zz',
              direc               = fs_out_dir,
              umin                = N_zz_min,
              umax                = N_zz_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
                                  
                                  
plot_variable(u                   = fs_srfmodel.M_ii,
              name                = 'M_ii',
              direc               = fs_out_dir,
              umin                = M_ii_min,
              umax                = M_ii_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = fs_srfmodel.M_ij,
              name                = 'M_ij',
              direc               = fs_out_dir,
              umin                = M_ij_min,
              umax                = M_ij_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = fs_srfmodel.M_iz,
              name                = 'M_iz',
              direc               = fs_out_dir,
              umin                = M_iz_min,
              umax                = M_iz_max,
              cmap                = n_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = fs_srfmodel.M_ji,
              name                = 'M_ji',
              direc               = fs_out_dir,
              umin                = M_ji_min,
              umax                = M_ji_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = fs_srfmodel.M_jj,
              name                = 'M_jj',
              direc               = fs_out_dir,
              umin                = M_jj_min,
              umax                = M_jj_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = fs_srfmodel.M_jz,
              name                = 'M_jz',
              direc               = fs_out_dir,
              umin                = M_jz_min,
              umax                = M_jz_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = fs_srfmodel.M_zi,
              name                = 'M_zi',
              direc               = fs_out_dir,
              umin                = M_zi_min,
              umax                = M_zi_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = fs_srfmodel.M_zj,
              name                = 'M_zj',
              direc               = fs_out_dir,
              umin                = M_zj_min,
              umax                = M_zj_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = fs_srfmodel.M_zz,
              name                = 'M_zz',
              direc               = fs_out_dir,
              umin                = M_zz_min,
              umax                = M_zz_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

#===============================================================================

plot_variable(u                   = rs_srfmodel.U3,
              name                = 'U_mag',
              umin                = U_min,
              umax                = U_max,
              vec_scale           = 22,
              direc               = rs_out_dir,
              cmap                = cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              cb_format           = '%g')

plot_variable(u                   = rs_bedmodel.p,
              name                = 'p',
              umin                = p_min,
              umax                = p_max,
              direc               = rs_out_dir,
              cmap                = cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              cb_format           = '%.1e')

plot_variable(u                   = rs_srfmodel.N_ii,
              name                = 'N_ii',
              direc               = rs_out_dir,
              umin                = N_ii_min,
              umax                = N_ii_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

plot_variable(u                   = rs_srfmodel.N_ij,
              name                = 'N_ij',
              direc               = rs_out_dir,
              umin                = N_ij_min,
              umax                = N_ij_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

plot_variable(u                   = rs_srfmodel.N_iz,
              name                = 'N_iz',
              direc               = rs_out_dir,
              umin                = N_iz_min,
              umax                = N_iz_max,
              cmap                = p_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

plot_variable(u                   = rs_srfmodel.N_ji,
              name                = 'N_ji',
              direc               = rs_out_dir,
              umin                = N_ji_min,
              umax                = N_ji_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

plot_variable(u                   = rs_srfmodel.N_jj,
              name                = 'N_jj',
              direc               = rs_out_dir,
              umin                = N_jj_min,
              umax                = N_jj_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

plot_variable(u                   = rs_srfmodel.N_jz,
              name                = 'N_jz',
              direc               = rs_out_dir,
              umin                = N_jz_min,
              umax                = N_jz_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = rs_srfmodel.N_zi,
              name                = 'N_zi',
              direc               = rs_out_dir,
              umin                = N_zi_min,
              umax                = N_zi_max,
              cmap                = p_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = rs_srfmodel.N_zj,
              name                = 'N_zj',
              direc               = rs_out_dir,
              umin                = N_zj_min,
              umax                = N_zj_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = rs_srfmodel.N_zz,
              name                = 'N_zz',
              direc               = rs_out_dir,
              umin                = N_zz_min,
              umax                = N_zz_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
                                  
                                  
plot_variable(u                   = rs_srfmodel.M_ii,
              name                = 'M_ii',
              direc               = rs_out_dir,
              umin                = M_ii_min,
              umax                = M_ii_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = rs_srfmodel.M_ij,
              name                = 'M_ij',
              direc               = rs_out_dir,
              umin                = M_ij_min,
              umax                = M_ij_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = rs_srfmodel.M_iz,
              name                = 'M_iz',
              direc               = rs_out_dir,
              umin                = M_iz_min,
              umax                = M_iz_max,
              cmap                = n_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = rs_srfmodel.M_ji,
              name                = 'M_ji',
              direc               = rs_out_dir,
              umin                = M_ji_min,
              umax                = M_ji_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = rs_srfmodel.M_jj,
              name                = 'M_jj',
              direc               = rs_out_dir,
              umin                = M_jj_min,
              umax                = M_jj_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = rs_srfmodel.M_jz,
              name                = 'M_jz',
              direc               = rs_out_dir,
              umin                = M_jz_min,
              umax                = M_jz_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = rs_srfmodel.M_zi,
              name                = 'M_zi',
              direc               = rs_out_dir,
              umin                = M_zi_min,
              umax                = M_zi_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = rs_srfmodel.M_zj,
              name                = 'M_zj',
              direc               = rs_out_dir,
              umin                = M_zj_min,
              umax                = M_zj_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = rs_srfmodel.M_zz,
              name                = 'M_zz',
              direc               = rs_out_dir,
              umin                = M_zz_min,
              umax                = M_zz_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

#===============================================================================

plot_variable(u                   = bp_srfmodel.U3,
              name                = 'U_mag',
              umin                = U_min,
              umax                = U_max,
              vec_scale           = 22,
              direc               = bp_out_dir,
              cmap                = cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              cb_format           = '%g')

plot_variable(u                   = bp_bedmodel.p,
              name                = 'p',
              umin                = p_min,
              umax                = p_max,
              direc               = bp_out_dir,
              cmap                = cmap,
              tp                  = True,
              show                = False,
              extend              = 'min',
              cb_format           = '%.1e')

plot_variable(u                   = bp_srfmodel.N_ii,
              name                = 'N_ii',
              direc               = bp_out_dir,
              umin                = N_ii_min,
              umax                = N_ii_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

plot_variable(u                   = bp_srfmodel.N_ij,
              name                = 'N_ij',
              direc               = bp_out_dir,
              umin                = N_ij_min,
              umax                = N_ij_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

plot_variable(u                   = bp_srfmodel.N_iz,
              name                = 'N_iz',
              direc               = bp_out_dir,
              umin                = N_iz_min,
              umax                = N_iz_max,
              cmap                = p_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

plot_variable(u                   = bp_srfmodel.N_ji,
              name                = 'N_ji',
              direc               = bp_out_dir,
              umin                = N_ji_min,
              umax                = N_ji_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

plot_variable(u                   = bp_srfmodel.N_jj,
              name                = 'N_jj',
              direc               = bp_out_dir,
              umin                = N_jj_min,
              umax                = N_jj_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')

plot_variable(u                   = bp_srfmodel.N_jz,
              name                = 'N_jz',
              direc               = bp_out_dir,
              umin                = N_jz_min,
              umax                = N_jz_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = bp_srfmodel.N_zi,
              name                = 'N_zi',
              direc               = bp_out_dir,
              umin                = N_zi_min,
              umax                = N_zi_max,
              cmap                = p_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = bp_srfmodel.N_zj,
              name                = 'N_zj',
              direc               = bp_out_dir,
              umin                = N_zj_min,
              umax                = N_zj_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = bp_srfmodel.N_zz,
              name                = 'N_zz',
              direc               = bp_out_dir,
              umin                = N_zz_min,
              umax                = N_zz_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
                                  
                                  
plot_variable(u                   = bp_srfmodel.M_ii,
              name                = 'M_ii',
              direc               = bp_out_dir,
              umin                = M_ii_min,
              umax                = M_ii_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = bp_srfmodel.M_ij,
              name                = 'M_ij',
              direc               = bp_out_dir,
              umin                = M_ij_min,
              umax                = M_ij_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = bp_srfmodel.M_iz,
              name                = 'M_iz',
              direc               = bp_out_dir,
              umin                = M_iz_min,
              umax                = M_iz_max,
              cmap                = n_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = bp_srfmodel.M_ji,
              name                = 'M_ji',
              direc               = bp_out_dir,
              umin                = M_ji_min,
              umax                = M_ji_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = bp_srfmodel.M_jj,
              name                = 'M_jj',
              direc               = bp_out_dir,
              umin                = M_jj_min,
              umax                = M_jj_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = bp_srfmodel.M_jz,
              name                = 'M_jz',
              direc               = bp_out_dir,
              umin                = M_jz_min,
              umax                = M_jz_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = bp_srfmodel.M_zi,
              name                = 'M_zi',
              direc               = bp_out_dir,
              umin                = M_zi_min,
              umax                = M_zi_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = bp_srfmodel.M_zj,
              name                = 'M_zj',
              direc               = bp_out_dir,
              umin                = M_zj_min,
              umax                = M_zj_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')
                                  
plot_variable(u                   = bp_srfmodel.M_zz,
              name                = 'M_zz',
              direc               = bp_out_dir,
              umin                = M_zz_min,
              umax                = M_zz_max,
              cmap                = sym_cmap,
              tp                  = True,
              show                = False,
              extend              = 'both',
              hide_ax_tick_labels = True,
              cb_format           = '%.1e')







