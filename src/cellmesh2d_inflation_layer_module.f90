!Cell Mesh 2D Inflation Layer Module 
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 0.1
!Updated 19-02-2024

!Module
module cellmesh2d_inflation_layer_mod
use io_utilities
use cellmesh2d_surface_mod
contains 


!Geometry inflation subroutine ================================================
subroutine inflate_geometry(surface_mesh,cm2dopt)
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt
type(surface_data) :: surface_mesh

!Variables - Local 
integer(in) :: ii!,jj
integer(in) :: n_inflayer
real(dp) :: infheight,cell_edgelength,nlevelR
real(dp) :: verticesC(surface_mesh%nvtx,2)
real(dp), dimension(:,:), allocatable :: grow_dir

!If the input number of layers is zero then guess appropriate number of layers from the highest refinement level cell size 
if (cm2dopt%inflayer_nlayer == 0) then 
    cell_edgelength = cm2dopt%far_field_bound/(2.0d0**(cm2dopt%Nrefine - 1))
    nlevelR = cm2dopt%inflayer_height/cell_edgelength
    cm2dopt%inflayer_nlayer = nint(nlevelR) + 1
    if (cm2dopt%dispt == 1) then
        write(*,'(A,I0,A)') '    {assinging ',cm2dopt%inflayer_nlayer,' layers within the inflation layer}'
    end if
end if 

!Extract base inflation layer height and number of layers
infheight = cm2dopt%inflayer_height 
n_inflayer = cm2dopt%inflayer_nlayer

!Store initial geometry
if (allocated(surface_mesh%vertices0)) then 
    deallocate(surface_mesh%vertices0)
end if 
allocate(surface_mesh%vertices0(surface_mesh%nvtx,2))
surface_mesh%vertices0(:,:) = surface_mesh%vertices(:,:)

!Allocate layer array 
if (allocated(surface_mesh%vertices_layer)) then 
    deallocate(surface_mesh%vertices_layer)
end if 
allocate(surface_mesh%vertices_layer(surface_mesh%nvtx,2,n_inflayer+1))

!Initialise growth direction array
allocate(grow_dir(surface_mesh%nvtx,2))
grow_dir(:,:) = 0.0d0 

!Inflate
verticesC(:,:) = surface_mesh%vertices(:,:)
do ii=1,n_inflayer

    !Display
    if (cm2dopt%dispt == 1) then
        write(*,'(A,I0)') '    level -> ', ii
    end if

    !Find growth directions to this layer
    call get_growth_direction(grow_dir,verticesC,surface_mesh,infheight/n_inflayer,cm2dopt)

    !Inflate to this layer with local step limit 
    verticesC(:,:) = verticesC(:,:) + grow_dir(:,:)

    !Store this layer
    surface_mesh%vertices_layer(:,:,n_inflayer-ii+1) = verticesC(:,:)

    !Debug -----
    ! open(11,file='io/vtxtestL'//int2str(ii))
    ! do jj=1,surface_mesh%nvtx
    !     write(11,*) verticesC(jj,:)
    ! end do 
    ! close(11)
end do 

!Store actual surface in the final layer
surface_mesh%vertices_layer(:,:,n_inflayer+1) = surface_mesh%vertices(:,:)

!Update surface to the highest layer
surface_mesh%vertices(:,:) = verticesC(:,:)
return
end subroutine inflate_geometry




!Build inflation layer mesh subroutine ================================================
subroutine mesh_inflation_layer(volume_mesh,surface_mesh,cm2dopt)
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ilayer
integer(in) :: Nsurfedge,Nsurfvtx,Ncell_new,Nedge_new,Nvertex_new,Nlayer
integer(in) :: surf_linkindex_temp(volume_mesh%nvtx)
integer(in), dimension(:), allocatable :: surface_edges,surface_vertices,surface_vertex_link,v2slist_vlink
integer(in), dimension(:), allocatable :: surfedge_ncell
integer(in), dimension(:,:), allocatable :: edge_new
real(dp), dimension(:), allocatable :: Cvol
real(dp), dimension(:,:), allocatable :: vertices_new

!Extract number of layers 
Nlayer = cm2dopt%inflayer_nlayer

!Find current surface mesh 
allocate(v2slist_vlink(volume_mesh%nvtx))
call get_surface_edges(Nsurfedge,surface_edges,volume_mesh%nedge,volume_mesh%edge)
call get_surface_vertices(Nsurfvtx,surface_vertices,v2slist_vlink,Nsurfedge,surface_edges,volume_mesh%nvtx,volume_mesh%edge)
deallocate(v2slist_vlink)

!Initialise new edges arrray 
Nedge_new = volume_mesh%nedge + (2*Nlayer + 1)*Nsurfedge
allocate(edge_new(Nedge_new,4))
edge_new(:,:) = 0 
edge_new(1:volume_mesh%nedge,:) = volume_mesh%edge(:,:)
Nedge_new = volume_mesh%nedge

!Initialise new vertex array 
Nvertex_new = volume_mesh%nvtx + (Nlayer + 1)*Nsurfvtx
allocate(vertices_new(Nvertex_new,2))
vertices_new(:,:) = 0 
vertices_new(1:volume_mesh%nvtx,:) = volume_mesh%vertices(:,:)

!Expand the linked surface index surf volume_mesh%surf_linkindex
surf_linkindex_temp(:) = volume_mesh%surf_linkindex(:)
deallocate(volume_mesh%surf_linkindex)
allocate(volume_mesh%surf_linkindex(Nvertex_new))
volume_mesh%surf_linkindex(:) = 0 
volume_mesh%surf_linkindex(1:volume_mesh%nvtx) = surf_linkindex_temp(:)

!Allocate volume to surface list linking array for vertices
allocate(v2slist_vlink(Nvertex_new))
v2slist_vlink(:) = 0 

!Reset vertex count 
Nvertex_new = volume_mesh%nvtx

!Initialise new cell count 
Ncell_new = volume_mesh%ncell

!Build each layer
do ilayer=1,Nlayer

    !Find current surface mesh 
    call get_surface_edges(Nsurfedge,surface_edges,Nedge_new,edge_new)
    call get_surface_vertices(Nsurfvtx,surface_vertices,v2slist_vlink,Nsurfedge,surface_edges,Nvertex_new,edge_new)

    !Add new vertex for each current surface vertex
    call add_layer_vertices(Nvertex_new,vertices_new,Nsurfvtx,surface_vertices,surface_vertex_link,&
    surface_mesh,volume_mesh,ilayer+1)

    !Build new cells for each surface edge
    call build_inflationonsurfedge_cells(surfedge_ncell,Ncell_new,Nsurfedge)

    !Build new layer internal edges
    call build_new_internal_edges(Nedge_new,edge_new,Nsurfvtx,surface_vertices,surface_vertex_link,&
    v2slist_vlink,Nsurfedge,surface_edges,surfedge_ncell)

    !Build new surface edges
    call build_new_surface_edges(Nedge_new,edge_new,surface_vertex_link,v2slist_vlink,Nsurfedge,surface_edges,surfedge_ncell) 

    !Update mesh surface link vertex indecies 
    call update_vslinks(volume_mesh,Nsurfvtx,surface_vertices,surface_vertex_link)

    !Display
    if (cm2dopt%dispt == 1) then
        write(*,'(A,I0)') '    level -> ', ilayer
    end if 
end do

!Update the actual volume mesh structure to include the inflation layer
volume_mesh%nvtx = Nvertex_new
deallocate(volume_mesh%vertices)
allocate(volume_mesh%vertices(volume_mesh%nvtx,2))
volume_mesh%vertices(:,:) = vertices_new(1:Nvertex_new,:)
volume_mesh%nedge = Nedge_new
deallocate(volume_mesh%edge)
allocate(volume_mesh%edge(volume_mesh%nedge,4))
volume_mesh%edge(:,:) = edge_new(1:Nedge_new,:)
volume_mesh%ncell = Ncell_new

!Evaluate volumes 
call get_cell_volumes(Cvol,volume_mesh)

!Completion of mesh construction display ----------------------
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '--> mesh construction completed'
    write(*,'(A,I0,A)') '    {cells: ',volume_mesh%ncell,'}'
    write(*,'(A,I0,A)') '    {edges: ',volume_mesh%nedge,'}'
    write(*,'(A,I0,A)') '    {vertices: ',volume_mesh%nvtx,'}'
    write(*,'(A,I0,A)') '    {surface intersection vertices: ',volume_mesh%nvtx_surf,'}'
    write(*,'(A,E11.5,A,E11.5,A,E11.5,A)') '    {cell volume (max/min/total) : ',maxval(Cvol),' / ',minval(Cvol),' / ',sum(Cvol),'}'
end if

!Debug -----
! open(11,file='io/vtxfull')
! do ilayer=volume_mesh%nvtx,Nvertex_new
!     write(11,*) vertices_new(ilayer,:)
! end do 
! close(11)
return 
end subroutine mesh_inflation_layer



!Build new internal edges subroutine ================================================
subroutine build_new_internal_edges(Nedge,edges,Nsurfvtx,surface_vertices,surface_vertex_link,&
                                    v2slist_vlink,Nsurfedge,surface_edges,surfedge_ncell)
implicit none 

!Variables - Import
integer(in) :: Nedge,Nsurfvtx,Nsurfedge
integer(in), dimension(:) :: surface_vertices,surface_edges,v2slist_vlink,surfedge_ncell
integer(in), dimension(:), allocatable :: surface_vertex_link
integer(in), dimension(:,:) :: edges

!Variables - Local 
integer(in) :: ii 
integer(in) :: v1,v2,c3,c4,etgt
integer(in) :: vtx_osedge_pn(Nsurfvtx,2)

!Find next and previous old surface edges for each old surface vertex
vtx_osedge_pn(:,:) = 0 !previous | next
do ii=1,Nsurfedge

    !Edge
    etgt = surface_edges(ii)

    !Vertices on this edge 
    v1 = edges(etgt,1)
    v2 = edges(etgt,2)

    !Indecies in the surface list 
    v1 = v2slist_vlink(v1)
    v2 = v2slist_vlink(v2)

    !Add this edge to each vertex
    vtx_osedge_pn(v1,2) = ii
    vtx_osedge_pn(v2,1) = ii
end do 

!Add new edge for each surface vertex 
do ii=1,Nsurfvtx

    !Vertices on this new edge
    v1 = surface_vertices(ii)
    v2 = surface_vertex_link(ii)

    !Find next and previous cells to link 
    if (vtx_osedge_pn(ii,1) .GT. 0) then 
        c3 = surfedge_ncell(vtx_osedge_pn(ii,1))
    else
        c3 = -2
    end if 
    if (vtx_osedge_pn(ii,2) .GT. 0) then 
        c4 = surfedge_ncell(vtx_osedge_pn(ii,2))
    else
        c4 = -2
    end if 

    !Build new edge
    Nedge = Nedge + 1
    edges(Nedge,1) = v1
    edges(Nedge,2) = v2

    !Link adjacent cells 
    edges(Nedge,3) = c3
    edges(Nedge,4) = c4
end do 
return 
end subroutine build_new_internal_edges




!Build new surface edges subroutine ================================================
subroutine build_new_surface_edges(Nedge,edges,surface_vertex_link,v2slist_vlink,Nsurfedge,surface_edges,surfedge_ncell)
implicit none 

!Variables - Import
integer(in) :: Nedge,Nsurfedge
integer(in), dimension(:) :: surface_edges,v2slist_vlink,surfedge_ncell
integer(in), dimension(:), allocatable :: surface_vertex_link
integer(in), dimension(:,:) :: edges

!Variables - Local 
integer(in) :: ii 
integer(in) :: edg0,v10,v20,v11,v21,ncell

!Construct a new edge for each existing surface edge
do ii=1,Nsurfedge

    !Existing surface edge to replace
    edg0 = surface_edges(ii)
    v10 = edges(edg0,1)
    v20 = edges(edg0,2)

    !Indecies of the vertices in the surface vertex list
    v10 = v2slist_vlink(v10)
    v20 = v2slist_vlink(v20)

    !New vertex indecies 
    v11 = surface_vertex_link(v10)
    v21 = surface_vertex_link(v20)

    !New cell to add
    ncell = surfedge_ncell(ii)

    !Build new edge
    Nedge = Nedge + 1
    edges(Nedge,1) = v11
    edges(Nedge,2) = v21

    !Set cell adjacency on new and existing edge
    if (edges(edg0,3) == -1) then 
        edges(Nedge,3) = -1
        edges(Nedge,4) = ncell 
        edges(edg0,3) = ncell 
    elseif (edges(edg0,4) == -1) then 
        edges(Nedge,3) = ncell 
        edges(Nedge,4) = -1
        edges(edg0,4) = ncell 
    else
        print *, '** warning: non surface edge detected within inflation layer growth'
    end if 
end do 
return 
end subroutine build_new_surface_edges




!Build new surface edge internal cells subroutine ================================================
subroutine build_inflationonsurfedge_cells(surfedge_ncell,Ncell,Nsurfedge)
implicit none 

!Variables - Import
integer(in) :: Ncell,Nsurfedge
integer(in), dimension(:), allocatable :: surfedge_ncell

!Variables - Local 
integer(in) :: ii 

!Reallocate link array
if (allocated(surfedge_ncell)) then 
    deallocate(surfedge_ncell)
end if 
allocate(surfedge_ncell(Nsurfedge))

!Assign cells
do ii=1,Nsurfedge
    Ncell = Ncell + 1
    surfedge_ncell(ii) = Ncell
end do 
return
end subroutine build_inflationonsurfedge_cells




!Update vertex surface links subroutine ================================================
subroutine update_vslinks(volume_mesh,Nsurfvtx,surface_vertices,surface_vertex_link)
implicit none 

!Variables - Import
integer(in) :: Nsurfvtx
integer(in), dimension(:) :: surface_vertices
integer(in), dimension(:), allocatable :: surface_vertex_link
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ii 
integer(in) :: vtgt0,vtgt1,lnkv

!Update the link for each surface vertex
do ii=1,Nsurfvtx

    !Original vertex index in the volume mesh
    vtgt0 = surface_vertices(ii)

    !New vertex index in the volume mesh 
    vtgt1 = surface_vertex_link(ii)

    !Update
    lnkv = volume_mesh%surf_linkindex(vtgt0)
    volume_mesh%surf_linkindex(vtgt0) = 0 
    volume_mesh%surf_linkindex(vtgt1) = lnkv
    !print *, lnkv,' || ',vtgt0,' -> ',vtgt1
end do 
return 
end subroutine update_vslinks




!Add vertices on layer subroutine ================================================
subroutine add_layer_vertices(Nvertex,vertices,Nsurfvtx,surface_vertices,surface_vertex_link,&
                              surface_mesh,volume_mesh,level)
implicit none 

!Variables - Import
integer(in) :: Nvertex,Nsurfvtx,level
integer(in), dimension(:) :: surface_vertices
integer(in), dimension(:), allocatable :: surface_vertex_link
real(dp), dimension(:,:) :: vertices 
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ii 
integer(in) :: vtgt,lidx,sidx,v1,v2
real(dp) :: sfrac
real(dp) :: vtxT(2)

!For each surface_vertices find its location through the surface link on the next level surface and insert this vertex
if (allocated(surface_vertex_link)) then 
    deallocate(surface_vertex_link)
end if 
allocate(surface_vertex_link(Nsurfvtx))
surface_vertex_link(:) = 0 
do ii=1,Nsurfvtx

    !Target vertex
    vtgt = surface_vertices(ii)

    !Linked surface index
    lidx = volume_mesh%surf_linkindex(vtgt)

    !Surface segment and fraction
    sidx = volume_mesh%surf_vtx_seg(lidx)
    sfrac = volume_mesh%surf_vtx_segfrac(lidx)
    v1 = surface_mesh%faces(sidx,1)
    v2 = surface_mesh%faces(sidx,2)

    !Vertex position 
    vtxT(:) = sfrac*surface_mesh%vertices_layer(v2,:,level) + (1.0d0 - sfrac)*surface_mesh%vertices_layer(v1,:,level)

    !Store and link this vertex
    Nvertex = Nvertex + 1
    vertices(Nvertex,:) = vtxT(:)
    surface_vertex_link(ii) = Nvertex
end do 
return 
end subroutine add_layer_vertices




!Subroutine to extract current surface mesh edges ================================================
subroutine get_surface_edges(Nsurfedge,surface_edges,Nedge,edges)
implicit none 

!Variables - Import
integer(in) :: Nsurfedge,Nedge
integer(in), dimension(:,:) :: edges 
integer(in), dimension(:), allocatable :: surface_edges

!Variables - Local 
integer(in) :: ii 

!Initialise 
Nsurfedge = 0 
if (allocated(surface_edges)) then 
    deallocate(surface_edges)
end if 

!Count and allocate surface edges 
Nsurfedge = 0 
do ii=1,Nedge
    if ((edges(ii,3) == -1) .OR. (edges(ii,4) == -1)) then 
        Nsurfedge = Nsurfedge + 1
    end if 
end do 
allocate(surface_edges(Nsurfedge))
surface_edges(:) = 0 
Nsurfedge = 0 

!Find surface edges 
do ii=1,Nedge
    if ((edges(ii,3) == -1) .OR. (edges(ii,4) == -1)) then 
        Nsurfedge = Nsurfedge + 1
        surface_edges(Nsurfedge) = ii 
    end if 
end do 
return 
end subroutine get_surface_edges




!Get surface vertices ================================================
subroutine get_surface_vertices(Nsurfvtx,surface_vertices,v2slist_vlink,Nsurfedge,surface_edges,Nvertex,edges)
implicit none 

!Variables - Import
integer(in) :: Nsurfedge,Nsurfvtx,Nvertex
integer(in), dimension(:), allocatable :: surface_edges,surface_vertices,v2slist_vlink
integer(in), dimension(:,:) :: edges 

!Variables - Local 
integer(in) :: ii 
integer(in) :: vtx_tag(Nvertex)

!Initialise 
Nsurfvtx = 0 
if (allocated(surface_vertices)) then 
    deallocate(surface_vertices)
end if 

!Count and allocate surface vertices 
vtx_tag(:) = 0 
do ii=1,Nsurfedge
    vtx_tag(edges(surface_edges(ii),1:2)) = 1
end do 
Nsurfvtx = sum(vtx_tag(:))
allocate(surface_vertices(Nsurfvtx))
surface_vertices(:) = 0 
Nsurfvtx = 0 

!Extract surface vertices 
v2slist_vlink(:) = 0
do ii=1,Nvertex
    if (vtx_tag(ii) == 1) then 
        Nsurfvtx = Nsurfvtx + 1
        surface_vertices(Nsurfvtx) = ii
        v2slist_vlink(ii) = Nsurfvtx
    end if 
end do 
return 
end subroutine get_surface_vertices




!Evaluate growth direction subroutine ================================================
subroutine get_growth_direction(grow_dir,verticesC,surface_mesh,steplen,cm2dopt)
implicit none 

!System data
type(cm2d_options) :: cm2dopt
type(surface_data) :: surface_mesh

!Variables - Import
real(dp) :: steplen
real(dp) ::  verticesC(surface_mesh%nvtx,2)
real(dp), dimension(:,:), allocatable :: grow_dir

!Variables - Local 
logical :: is_inverted
integer(in) :: ii,ff,vv,ss
integer(in) :: nv,pv,max_int_iter,v1,v2
real(dp) :: sstep,dstep,stepstab,normv1,normv2,dnorm,gweight,alpha,h_total,h_total0,h_totalP,h_stepav
real(dp) :: verticesP(surface_mesh%nvtx,2),verticesT(surface_mesh%nvtx,2)
real(dp) :: df(surface_mesh%nfcs),dv(surface_mesh%nvtx),stepnorm(surface_mesh%nvtx),vec1(2),vec2(2)
real(dp), dimension(:,:), allocatable :: grow_dirL

!Set maximum intergration iterations
max_int_iter = cm2dopt%inflayer_nintstep_max

!Set evening stepsize
sstep = cm2dopt%inflayer_eveningstepsize

!Set dv step perturbation 
dstep = cm2dopt%inflayer_dvstepsize

!Set stabilisation parameter
stepstab = 1e-8

!Initialise as normals
call evaluate_surface_normals(grow_dirL,surface_mesh,verticesC)
if (cm2dopt%meshinout == 'in') then 
    grow_dirL(:,:) = -grow_dirL(:,:) 
end if 

!Integrate along step
h_total = 0.0d0 
verticesP(:,:) = verticesC(:,:)
do ii=1,max_int_iter

    !Evaluate df and dv for each face 
    dv(:) = 0.0d0 
    do ff=1,surface_mesh%nfcs
        v1 = surface_mesh%faces(ff,1)
        v2 = surface_mesh%faces(ff,2)
        vec1(:) = verticesP(v2,:) - verticesP(v1,:)
        vec2(:) = verticesP(v1,:) - verticesP(v2,:)
        df(ff) = dot_product(grow_dirL(v1,:),vec1(:)) + dot_product(grow_dirL(v2,:),vec2(:))
        dv(v1) = dv(v1) + df(ff)
        dv(v2) = dv(v2) + df(ff)
    end do 
    dv(:) = 0.5d0*dv(:)
    do vv=1,surface_mesh%nvtx
        ! dv(vv) = dv(vv)/(1.0d0 + exp(-500.0d0*dv(vv)))
        if (dv(vv) .LT. 0.0d0) then 
            ! dv(vv) = cm2dopt%inflayer_cvxep*dv(vv)
            dv(vv) = cm2dopt%inflayer_cvxdp*dv(vv)
        end if 
    end do 

    !Step vertices to valid configuration with no crossovers 
    alpha = 1.0d0 
    do ss=1,cm2dopt%inflayer_nbclinesearch

        !Test step
        do vv=1,surface_mesh%nvtx

            !Find growth weighting 
            gweight = (1.0d0 + dstep*dv(vv))
            ! gweight = max(min(gweight/(abs(surface_mesh%vtx_rcurv(vv)) + stepstab),1.0d0),gweight)
            if (gweight .LT. 0.0d0) then 
                gweight = 0.0d0 
            end if 

            !Add damping near sharp vertices
            if (surface_mesh%vtx_sharp(vv) == 1) then 
                gweight = 1.0d0
            end if 

            !Test step vertices 
            verticesT(vv,1) = verticesP(vv,1) + alpha*steplen*grow_dirL(vv,1)*gweight
            verticesT(vv,2) = verticesP(vv,2) + alpha*steplen*grow_dirL(vv,2)*gweight
        end do 

        !Check surface 
        is_inverted = has_face_inversion(verticesT,verticesP,surface_mesh)
        if (is_inverted) then !if inverted somewhere then reduct stepsize 
            alpha = 0.5d0*alpha
        else !Accept step if valid surface
            exit 
        end if 
    end do 

    !Find average step length h_total
    stepnorm(:) = sqrt((verticesT(:,1) - verticesP(:,1))**2 + (verticesT(:,2) - verticesP(:,2))**2)
    h_stepav = sum(stepnorm(:))/real(surface_mesh%nvtx,dp)

    !Predict total layer height 
    h_total0 = h_total
    h_totalP = h_total + h_stepav

    !Damp step to maximum layer height if required
    alpha = 1.0d0 
    if (h_totalP .GE. steplen) then 
        alpha = (steplen - h_total0)/h_stepav
    end if 

    !Set new vertex positions
    verticesP(:,:) = verticesP(:,:) + alpha*(verticesT(:,:) - verticesP(:,:) )

    !Accumulate total layer height
    h_total = h_total0 + alpha*h_stepav

    !Evening smooth along surface
    do vv=1,surface_mesh%nvtx
        if (surface_mesh%vtx_sharp(vv) == 0) then 

            !Next and previous vertices
            pv = get_previous_vertex(surface_mesh,vv)
            nv = get_next_vertex(surface_mesh,vv)

            !Adjacent vertex vectors
            vec1(:) = verticesP(pv,:) - verticesP(vv,:)
            vec2(:) = verticesP(nv,:) - verticesP(vv,:)

            !Evening norm
            normv1 = norm2(vec1)
            normv2 = norm2(vec2)
            dnorm = abs(normv2 - normv1)/(normv1 + normv2 + stepstab)
            if (dv(vv) .LT. 0.0d0) then 
                dnorm = cm2dopt%inflayer_cvxep*dnorm
            end if 

            !Even vertex
            ! if (normv1 .GT. normv2) then 
            !     verticesT(vv,:) = verticesP(vv,:) + vec1(:)*dnorm*sstep
            ! elseif (normv2 .GT. normv1) then 
            !     verticesT(vv,:) = verticesP(vv,:) + vec2(:)*dnorm*sstep
            ! else
            !     verticesT(vv,:) = verticesP(vv,:)
            ! end if 
            verticesT(vv,:) = verticesP(vv,:) + 0.5d0*dnorm*sstep*(vec1(:) + vec2(:))
        else
            verticesT(vv,:) = verticesP(vv,:)
        end if 
    end do 
    verticesP(:,:) = verticesT(:,:)

    !Update local direction
    call evaluate_surface_normals(grow_dirL,surface_mesh,verticesP)
    if (cm2dopt%meshinout == 'in') then 
        grow_dirL(:,:) = -grow_dirL(:,:) 
    end if 

    !Exit if approximate layer height has been reached by htotal 
    if (h_total .GE. steplen) then 
        if (cm2dopt%dispt == 1) then
            write(*,'(A,I0,A,A,A)') '    {complete at niter = ',ii,' with height = ',real2F0_Xstring(h_total,8_in),'}'
        end if 
        exit 
    end if 

    !Debug -----
    ! open(11,file='io/vtxtest'//int2str(ii))
    ! do ff=1,surface_mesh%nvtx
    !     write(11,*) verticesP(ff,:)
    ! end do 
    ! close(11)
    ! open(11,file='io/growdir'//int2str(ii))
    ! do ff=1,surface_mesh%nvtx
    !     write(11,*) grow_dirL(ff,:)
    ! end do 
    ! close(11)
end do 

!Evaluate overall direction 
grow_dir(:,:) = verticesP(:,:) - verticesC(:,:)

! open(11,file='io/vtx_level')
! do ff=1,surface_mesh%nvtx
!     write(11,*) verticesP(:,:)
! end do 
! close(11)
! open(11,file='io/vtx_level')
! do ff=1,surface_mesh%nvtx
!     write(11,*) verticesP(ff,:) !+ grow_dir(ff,:)
! end do 
! close(11)

return 
end subroutine get_growth_direction




!Function to check for face inversion anywhere on the surface ================================================
function has_face_inversion(verticesN,vertices0,surface_mesh) result(invertstate)
implicit none 

!Result
logical :: invertstate

!Variables - Import
real(dp), dimension(:,:) :: verticesN,vertices0
type(surface_data) :: surface_mesh

!Variables - Local 
integer(in) :: ff 
integer(in) :: v1,v2
real(dp) :: face0(2),faceN(2)

!Set initial state
invertstate = .false.

!Check for inverted faces
do ff=1,surface_mesh%nfcs

    !Vertices
    v1 = surface_mesh%faces(ff,1)
    v2 = surface_mesh%faces(ff,2)

    !Original and new face vectors
    face0(:) = vertices0(v2,:) - vertices0(v1,:)
    faceN(:) = verticesN(v2,:) - verticesN(v1,:)

    ! print *, dot_product(face0,faceN)

    !Check for inversion
    if (dot_product(face0,faceN) .LE. 0.0d0) then 
        invertstate = .true.
        exit 
    end if 
end do 
return 
end function has_face_inversion


end module cellmesh2d_inflation_layer_mod