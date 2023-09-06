!cell_mesh2d data types module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 0.6
!Updated 15-08-2023

!Module
module cellmesh2d_data_mod

!Integer data types 
use ISO_FORTRAN_ENV, only: in32=>int32
use ISO_FORTRAN_ENV, only: in=>int64

!Double precision data type 
use ISO_FORTRAN_ENV, only: dp=>real64     

!Set round off error bound 
real(dp), parameter :: min_precision = 1E-12_dp 

!Options data type
type cm2d_options
    character(len=:), allocatable :: iopath,optpath,surfacename,surface_dir,boundary_dir,meshinout,meshfrmat
    integer(in) :: Nrefine,NrefineB,Ncell_max,Nrefine_flood_i,Nrefine_flood_f,Nrefine_flood_b,meshtype
    integer(in) :: dispt,NintEmax,glink_con,glink_nnn,glink_nsmooth,ADTmax_depth,nlpflood,nlpsmooth,set_mbounds
    integer(in) :: set_customBCs,remFFzones,remISzones,remNCzones,Nsstype,Nzone_cBC,NPsinterp,surface_type,cm2dfailure
    real(dp) :: ADTpadding,far_field_bound,EminLength,CminVol,elenpad,intcointol,surfRcurvM,om_offset_x,om_offset_y
    real(dp) :: mxmin,mxmax,mymin,mymax
    integer(in), dimension(:), allocatable :: BC_zone_bc
    real(dp), dimension(:,:), allocatable :: BC_zone_coords
end type cm2d_options

!Nolume mesh data type 
type vol_mesh_data
    integer(in) :: nvtx,nedge,ncell,nvtx_surf
    integer(in), dimension(:), allocatable :: cell_level,cell_qtidx,vtx_type,vtx_surfseg
    integer(in), dimension(:), allocatable :: surf_vtx,surf_vtx_seg
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
    integer(in), dimension(:,:), allocatable :: faces,v2f,smvtx2vmedge
    real(dp), dimension(:), allocatable :: face_rcurv,vtx_rcurv,vtx_rsearch
    real(dp), dimension(:,:), allocatable :: vertices,vertices_full
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
    integer(in), dimension(:), allocatable :: vertices,edges
end type cell_data

end module cellmesh2d_data_mod