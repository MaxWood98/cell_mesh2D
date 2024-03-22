!cell_mesh2d data types module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 1.0
!Updated 22-03-2024

!Module
module cellmesh2d_data_mod

!Integer data types 
use ISO_FORTRAN_ENV, only: in32=>int32
use ISO_FORTRAN_ENV, only: in=>int64

!Double precision data type 
use ISO_FORTRAN_ENV, only: dp=>real64     

!Set round off error bound 
real(dp), parameter :: min_precision = 1E-12_dp 

!Surface from volume interpolation type
type surfFvolinterp
    integer(in) :: npnts_vi
    integer(in), dimension(:), allocatable :: vol_pnts,surf_smooth_pnts
    real(dp), dimension(:), allocatable :: surf2volRBF,surf_smoothRBF
    real(dp), dimension(:,:), allocatable :: Ri
end type surfFvolinterp

!Cell data type 
type cell_data 
    integer(in) :: nvtx,nedge 
    integer(in), dimension(:), allocatable :: vertices,edges,boundary_condition
end type cell_data

!Options data type
type cm2d_options
    logical :: dispt
    character(len=:), allocatable :: mode,iopath,optpath,surface_filename,options_filename,bcondzone_filename
    character(len=:), allocatable :: surface_dir,boundary_dir,meshinout,meshtype,surface_type,Nsstype,glink_con_exp
    character(len=:), allocatable :: glink_type,build_inflayer,srfinclean,set_mbounds,set_customBCs,remFFzones,remISzones,remNCzones
    integer(in) :: Nrefine,NrefineB,Ncell_max,Nrefine_flood_i,Nrefine_flood_f,Nrefine_flood_b,inflayer_nlayer
    integer(in) :: inflayer_nintstep_max,inflayer_nbclinesearch,dfield_niter,inflayer_enflood,inflayer_ensubiter
    integer(in) :: inflayer_stepnacheck
    integer(in) :: NintEmax,glink_nnn,glink_nsmooth,ADTmax_depth,ADTminNodedivsize,nlpflood,nlpsmooth
    integer(in) :: Nzone_cBC,NPsinterp,cm2dfailure
    real(dp) :: ADTpadding,far_field_bound,sEminLength,EminLength,CminVol,elenpad,intcointol,surfRcurvM,om_offset_x,om_offset_y
    real(dp) :: mxmin,mxmax,mymin,mymax,RBF_relaxP,RBF_relaxD,RBF_rsup,vtx_sharp_dpval
    real(dp) :: dfield_cfl,dfield_kd,dfield_cres,inflayer_enormw
    real(dp) :: inflayer_height,inflayer_ew,inflayer_dvw,inflayer_h0,inflayer_cvxep,inflayer_cvxdp,inflayer_lreb,inflayer_ebcbase
    integer(in), dimension(:), allocatable :: BC_zone_bc
    real(dp), dimension(:,:), allocatable :: BC_zone_coords
end type cm2d_options

!Nolume mesh data type 
type vol_mesh_data
    integer(in) :: nvtx,nedge,ncell,nvtx_surf
    integer(in), dimension(:), allocatable :: cell_level,cell_qtidx,vtx_type,vtx_surfseg
    integer(in), dimension(:), allocatable :: surf_vtx,surf_vtx_seg,surf_linkindex
    integer(in), dimension(:,:), allocatable :: V2E
    integer(in), dimension(:,:), allocatable :: edge !v1 v2 cellL cellR
    real(dp), dimension(:), allocatable :: surf_vtx_segfrac
    real(dp), dimension(:,:), allocatable :: vertices,edge_midpoint
    type(surfFvolinterp), dimension(:), allocatable :: vsinterp
    type(cell_data), dimension(:), allocatable :: cells 
end type vol_mesh_data

!Surface data type 
type surface_data
    integer(in) :: nvtx,nfcs,nvtxf
    integer(in), dimension(:), allocatable :: vtx_active,vtx_sharp
    integer(in), dimension(:,:), allocatable :: faces,v2f,smvtx2vmedge
    real(dp), dimension(:), allocatable :: face_rcurv,vtx_rcurv,vtx_rsearch
    real(dp), dimension(:,:), allocatable :: vertices,vertices_full,normals,vertices0
    real(dp), dimension(:,:,:), allocatable :: vertices_layer
end type surface_data

!Quadtree data type 
type quadtree_data
    integer(in) :: cins,vins 
    integer(in), dimension(:), allocatable :: cell_level
    integer(in), dimension(:,:), allocatable :: cell_vcnr,cell_vemid,cell_child,cell_adjacent
    real(dp), dimension(:,:), allocatable :: vtx
end type quadtree_data

!Edge intersection data type
type edgeint 
    integer(in) :: type
    integer(in) :: nint,refend1
    integer(in), dimension(:), allocatable :: vtx_idx,inttype,intidx,surfseg,int_seg_type,int_inout
    integer(in), dimension(:), allocatable :: vncell,nfint,vmedge,vmeshcell
    integer(in), dimension(:,:), allocatable :: edge_mesh,vcell
    real(dp), dimension(:), allocatable :: intfrac
    real(dp), dimension(:,:), allocatable :: intloc 
end type edgeint 

!Mesh data data type 
type meshdata 
    integer(in), dimension(:), allocatable :: vtx_vlnc,cell_nedge,bc_active
    integer(in), dimension(:,:), allocatable :: vtx_2_cell,vtx_v2v,vtx_v2e,cell2edge
    real(dp), dimension(:), allocatable :: cell_area,edgedx,edgedy,edgedxd,edgedyd
    real(dp), dimension(:,:), allocatable :: cell_midp,edge_midp,cell_vtx_W
end type meshdata

end module cellmesh2d_data_mod