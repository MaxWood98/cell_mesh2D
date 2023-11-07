!Cell Mesh 2D Cutcell Mesh Generation Module 
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 1.8
!Updated 18-09-2023

!Module
module cellmesh2d_mesh_generation_mod
use cellmesh2d_io_mod
use cellmesh2d_surface_mod
use cellmesh2d_quadtree_mod
use cellmesh2d_postprocess_mod
use cellmesh2d_mesh_build_exact_mod
use cellmesh2d_gradient_coupling_mod
contains

!cell_mesh2d mesh construction subroutine ================================================
subroutine cell_mesh2d_mesh(volume_mesh,surface_mesh,cm2dopt)
implicit none 

!Variables - Import ----
type(cm2d_options) :: cm2dopt
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh

!Variables - Local ----
!Counters 
integer(in) :: ii

!Object data
real(dp) :: obj_max_x,obj_max_y,obj_min_x,obj_min_y,obj_cx,obj_cy
real(dp), dimension(:,:), allocatable :: tvtx

!Mesh data
integer(in) :: MaxValence,Nmerge,Nmerge_fail,nmiter,nifail
real(dp), dimension(:), allocatable :: Cvol

!AD_Tree data
integer(in) :: Ndim,node_minDIVsize
real(dp) :: global_target_pad
type(tree_data) :: surface_adtree

!Quadtree 
integer(in), dimension(:), allocatable :: cell_keep,vtx_external,vtx_type1
type(quadtree_data) :: qt_mesh

!Set hardcoded parameters -------------------------------------------------------
!Set MaxValence parameter (=4)
MaxValence = 4

!Set global object bounding box padding for adtree node containement
global_target_pad = 0.0d0

!Adtree number of dimensions (4)
Ndim = 4

!Set minimum divisible node size within the AD tree 
node_minDIVsize = 10

!Initialise failure tag
cm2dopt%cm2dfailure = 0

!Set initial count state
volume_mesh%nvtx = 0 
volume_mesh%nedge = 0 
volume_mesh%ncell = 0 

!Property calculation -------------------------------------------------------
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '--> constructing surface geometry parameters'
end if

!Preprocess surface mesh 
call preprocess_surface_mesh(surface_mesh,cm2dopt)
if (cm2dopt%cm2dfailure == 1) then
    return 
end if 

!Global object bounds
obj_max_x = maxval(surface_mesh%vertices(:,1))
obj_max_y = maxval(surface_mesh%vertices(:,2))
obj_min_x = minval(surface_mesh%vertices(:,1))
obj_min_y = minval(surface_mesh%vertices(:,2))

!Global object centre
obj_cx = 0.5d0*(obj_max_x + obj_min_x)
obj_cy = 0.5d0*(obj_max_y + obj_min_y)

!Displays
if (cm2dopt%dispt == 1) then
    write(*,'(A)') ' '
    write(*,'(A)') '== Properties =========================='
    write(*,"(A,I0)") '   Imported surface vertices = ', surface_mesh%nvtx
    write(*,"(A,F12.6,A,F12.6)") '   Geometry dimensions (dx/dy) = ', obj_max_x - obj_min_x,' ', obj_max_y - obj_min_y
    write(*,'(A)') '== Properties =========================='
    write(*,'(A)') ' '
end if
if (cm2dopt%cm2dfailure == 1) then
    return 
end if 

!Construct AD_Tree on target mesh ------------------------------------------
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '--> constructing AD-tree on target geometry surface mesh'
end if

!Construct 4D bounding box coordinates for each face in the surface mesh
allocate(tvtx(surface_mesh%nfcs,4))
do ii=1,surface_mesh%nfcs
    tvtx(ii,1) = minval(surface_mesh%vertices(surface_mesh%faces(ii,:),1)) !xmin
    tvtx(ii,2) = minval(surface_mesh%vertices(surface_mesh%faces(ii,:),2)) !ymin
    tvtx(ii,3) = maxval(surface_mesh%vertices(surface_mesh%faces(ii,:),1)) !xmax
    tvtx(ii,4) = maxval(surface_mesh%vertices(surface_mesh%faces(ii,:),2)) !ymax
end do

!Construct ad_tree
call build_ADtree(surface_adtree,ndim,cm2dopt%ADTmax_depth,node_minDIVsize,tvtx,global_target_pad,cm2dopt%dispt)

!Quadtree construction -------------------------------------------------------
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '--> constructing mesh quadtree'
end if

!Quadree refinement to geometry
call quadtree_mesh_refine(qt_mesh,surface_mesh,surface_adtree,cm2dopt,obj_cx,obj_cy)
if (cm2dopt%cm2dfailure == 1) then
    return 
end if 

!Identify geometry internal, external and overlapping cells
call quadtree_mesh_trim2geom(cell_keep,vtx_external,vtx_type1,qt_mesh,surface_mesh,surface_adtree,cm2dopt)
if (cm2dopt%cm2dfailure == 1) then
    return 
end if 

!Construct volume mesh -------------------------------------------------------
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '--> constructing full volume mesh'
end if
call construct_mesh_EXACT(volume_mesh,surface_mesh,cell_keep,vtx_external,vtx_type1,qt_mesh,surface_adtree,cm2dopt)
if (cm2dopt%cm2dfailure == 1) then 
    return 
end if

!Remove any double boundary condition edges from surfaces outside the mesh domain 
call remove_double_bc_edges(volume_mesh,cm2dopt)

!Return if no mesh edges constructed
if (volume_mesh%nedge == 0) then !Check for complete mesh removal 
    cm2dopt%cm2dfailure = 1
    print *, '** no mesh edges constructed after clipping to geometry'
    return 
end if

!Set sequential cell indecies 
call remap_cell_indecies(volume_mesh)

!Postprocess complete mesh ---------------------------------------------------
!Apply custom boundary conditions in target regions 
if (cm2dopt%set_customBCs == 1) then 
    call set_custom_bcs(volume_mesh,cm2dopt) 
end if 

!Remove mesh zones with any far field boundary conditions 
if (cm2dopt%remFFzones == 1) then 
    call remove_farfield_adjacent(volume_mesh) 
    if (volume_mesh%nedge == 0) then !Check for complete mesh removal 
        cm2dopt%cm2dfailure = 1
        print *, '** mesh partitioning error -> no retained edges [remove far field zones]'
        return 
    end if
end if

!Remove mesh zones with no custom boundary conditions 
if (cm2dopt%remNCzones == 1) then 
    call remove_ncb_adjacent(volume_mesh,cm2dopt)
    if (volume_mesh%nedge == 0) then !Check for complete mesh removal 
        cm2dopt%cm2dfailure = 1
        print *, '** mesh partitioning error -> no retained edges [remove not custom boundary condition connected zones]'
        return 
    end if
end if

!Remove mesh zones that are only connected to wall boundary conditions (isolated regions of the mesh with no inflow/outflow possible)
if (cm2dopt%remISzones == 1) then 
    call remove_isolated_regions(volume_mesh)
    if (volume_mesh%nedge == 0) then !Check for complete mesh removal 
        cm2dopt%cm2dfailure = 1
        print *, '** mesh partitioning error -> no retained edges [remove isolated zones]'
        return 
    end if
end if 

!Check for and correct bisected cells 
call correct_bisected_cells(volume_mesh,cm2dopt)

!Clean mesh by removing short edges and sliver cells 
call clean_mesh_shortE(volume_mesh,cm2dopt)
nmiter = 10 !volume_mesh%ncell !set maximum merging iterations 
nifail = 0 
do ii=1,nmiter
    call clean_mesh_sliverC(volume_mesh,cm2dopt,Nmerge,Nmerge_fail)
    if (Nmerge_fail .NE. 0) then 
        nifail = nifail + 1
    end if 
    if (Nmerge .NE. 0) then 
        nifail = nifail - 1
    end if 
    if ((Nmerge == 0) .AND. (Nmerge_fail .NE. 0) .AND. (nifail == 1)) then 
        exit 
    end if 
    if ((Nmerge == 0) .AND. (Nmerge_fail == 0)) then 
        exit 
    end if
end do 

!Remap cell indecies 
call remap_cell_indecies(volume_mesh)

!Remove mesh internal valence two vertices 
call clean_internal_vlnc2_vertices(volume_mesh,cm2dopt)

!Simplify the mesh surface within each cell if simplified surface requested
if (cm2dopt%surface_type == 0) then     
    if (cm2dopt%dispt == 1) then
        write(*,'(A)') '--> constructing simplified surface geometry'
    end if
    call simplify_surface(volume_mesh)
end if 

!Build surface vertex mappings 
call build_surface_links(volume_mesh,surface_mesh)

!Flip surface and boundary normals if requested (both are constructed in state 'in')
if (cm2dopt%surface_dir == 'out') then 
    call flip_set_edges(volume_mesh,-1_in)
end if 
if (cm2dopt%boundary_dir == 'out') then 
    call flip_set_edges(volume_mesh,-2_in)
    call flip_set_edges(volume_mesh,-3_in)
    call flip_set_edges(volume_mesh,-4_in)
    call flip_set_edges(volume_mesh,-5_in)
    call flip_set_edges(volume_mesh,-6_in)
end if 

!Restructure mesh for SU2 output formats 
if (cm2dopt%meshfrmat == 'su2_cutcell') then 
    call split_vlnc_gt4cells(volume_mesh)
    call remap_cell_indecies(volume_mesh)
    call build_mesh_cells(volume_mesh)
elseif (cm2dopt%meshfrmat == 'su2_dual') then 
    call construct_dual_mesh(volume_mesh)
    call deform_mesh2surface(volume_mesh,surface_mesh,surface_adtree)
    call remap_cell_indecies(volume_mesh)
    call build_mesh_cells(volume_mesh)
end if 

!Apply near surface mesh smoothing 
if ((cm2dopt%Nsstype == 1) .OR. (cm2dopt%meshtype == 1)) then 
    call nearsurf_lap_smooth(volume_mesh,cm2dopt,MaxValence)
end if 

!Evaluate cell volumes to check for invalid cells 
call get_cell_volumes(Cvol,volume_mesh)

!Construct gradient mesh to surface coupling matrix if requeted 
if (cm2dopt%glink_con == 1) then 
    if (cm2dopt%dispt == 1) then
        write(*,'(A)') '--> constructing volume-surface gradient coupling matrix'
    end if
    call construct_surfvol_grad_coupling(volume_mesh,surface_mesh,Ndim,node_minDIVsize,global_target_pad,cm2dopt)
end if 

!Completion of mesh construction display ----------------------
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '--> mesh construction completed'
    write(*,'(A,I0,A)') '    {cells: ',volume_mesh%ncell,'}'
    write(*,'(A,I0,A)') '    {edges: ',volume_mesh%nedge,'}'
    write(*,'(A,I0,A)') '    {vertices: ',volume_mesh%nvtx,'}'
    write(*,'(A,I0,A)') '    {surface intersection vertices: ',volume_mesh%nvtx_surf,'}'
    write(*,'(A,E11.5,A,E11.5,A,E11.5,A)') '    {cell volume (max/min/total) : ',maxval(Cvol),' / ',minval(Cvol),' / ',sum(Cvol),'}'
end if
if (minval(Cvol) .LE. 0.0d0) then 
    ! cm2dopt%cm2dfailure = 1
    do ii=1,volume_mesh%ncell
        if (Cvol(ii) .LT. 0.0d0) then 
            print '(A,I0)', '** negative volume cell identified: ',ii
        elseif (Cvol(ii) == 0.0d0) then  
            print '(A,I0)', '** zero volume cell identified: ',ii
        end if 
    end do 
end if  
return 
end subroutine cell_mesh2d_mesh




!Cell volumes subroutine ===========================
subroutine get_cell_volumes(Cvol,volume_mesh)
implicit none 

!Variables - Import
real(dp), dimension(:), allocatable :: Cvol
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ii
real(dp) :: AEdge

if (allocated(Cvol)) then 
    deallocate(Cvol)
end if
allocate(Cvol(volume_mesh%ncell))
Cvol(:) = 0.0d0
do ii=1,volume_mesh%nedge
    AEdge = Asegment(volume_mesh%vertices(volume_mesh%edge(ii,1),:),volume_mesh%vertices(volume_mesh%edge(ii,2),:))
    if (volume_mesh%edge(ii,4) .GT. 0) then
        Cvol(volume_mesh%edge(ii,4)) = Cvol(volume_mesh%edge(ii,4)) - AEdge
    end if 
    if (volume_mesh%edge(ii,3) .GT. 0) then
        Cvol(volume_mesh%edge(ii,3)) = Cvol(volume_mesh%edge(ii,3)) + AEdge
    end if
end do
return 
end subroutine get_cell_volumes


end module cellmesh2d_mesh_generation_mod