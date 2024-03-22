!Cell Mesh 2D Inflation Layer Module 
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 0.5
!Updated 23-02-2024

!Module
module cellmesh2d_inflation_layer_mod
use io_utilities
use cellmesh2d_surface_mod
use cellmesh2d_mesh_generation_mod
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
real(dp) :: even_bias(surface_mesh%nvtx)
real(dp) :: verticesC(surface_mesh%nvtx,2)
real(dp), dimension(:,:), allocatable :: grow_dir
type(tree_data) :: surface_adtree

!Construct adtree on base surface geometry 
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '--> constructing AD-tree on target geometry surface mesh'
end if
allocate(surface_adtree%tvtx(surface_mesh%nfcs,4))
do ii=1,surface_mesh%nfcs
    surface_adtree%tvtx(ii,1) = minval(surface_mesh%vertices(surface_mesh%faces(ii,:),1)) !xmin
    surface_adtree%tvtx(ii,2) = minval(surface_mesh%vertices(surface_mesh%faces(ii,:),2)) !ymin
    surface_adtree%tvtx(ii,3) = maxval(surface_mesh%vertices(surface_mesh%faces(ii,:),1)) !xmax
    surface_adtree%tvtx(ii,4) = maxval(surface_mesh%vertices(surface_mesh%faces(ii,:),2)) !ymax
end do

!Construct ad_tree
call build_ADtree(surface_adtree,4_in,cm2dopt%ADTmax_depth,cm2dopt%ADTminNodedivsize,surface_adtree%tvtx,0.0d0,cm2dopt%dispt)

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

!Initialise 
verticesC(:,:) = surface_mesh%vertices(:,:)

!Build evening bias array 
call construct_even_bias(even_bias,verticesC,surface_mesh,surface_adtree,cm2dopt)

!Inflate
do ii=1,n_inflayer

    !Display
    if (cm2dopt%dispt == 1) then
        write(*,'(A,I0)') '    level -> ', ii
    end if

    !Find growth directions to this layer
    call get_growth_direction(grow_dir,verticesC,even_bias,surface_mesh,infheight/n_inflayer,cm2dopt)

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




!Construct even bias subroutine ================================================
subroutine construct_even_bias(even_bias,verticesC,surface_mesh,surface_adtree,cm2dopt)
implicit none 

!System data
type(cm2d_options) :: cm2dopt
type(surface_data) :: surface_mesh
type(tree_data) :: surface_adtree

!Variables - Import
real(dp) :: even_bias(surface_mesh%nvtx),verticesC(surface_mesh%nvtx,2)

!Variables - Local 
integer(in) :: vv,nn,kk,ff
integer(in) :: nselected,ftgt,int_type,fx1,fx2,vp,vn
integer(in) :: node_select(surface_adtree%nnode),ebias_fixed(surface_mesh%nvtx)
real(dp) :: zxmin,zxmax,zymin,zymax,dv,stabval
real(dp) :: vtx1(2),vtx2(2),vtx1f(2),vtx2f(2),vecp(2),vecn(2)
real(dp) :: even_biasF(surface_mesh%nvtx)
real(dp), dimension(:,:), allocatable :: normals
    
!Set stabilisation value 
stabval = 1e-8

!Construct local normals 
call evaluate_surface_normals(normals,surface_mesh,verticesC)
if (cm2dopt%meshinout == 'in') then 
    normals(:,:) = -normals(:,:) 
end if 

!Initialise bias 
even_bias(:) = 0.0d0 

!Identify vertices for which evening will be allowed 
do vv=1,surface_mesh%nvtx !Add from normal intersection

    !Direction vector 
    vtx1(:) = verticesC(vv,:)
    vtx2(:) = verticesC(vv,:) + normals(vv,:)

    !Excluded faces
    fx1 = surface_mesh%v2f(vv,1)
    fx2 = surface_mesh%v2f(vv,2)

    !Check for intersection candidates
    zxmin = min(vtx1(1),vtx2(1))
    zxmax = max(vtx1(1),vtx2(1))
    zymin = min(vtx1(2),vtx2(2))
    zymax = max(vtx1(2),vtx2(2))
    call search_ADtree(nselected,node_select,surface_adtree,zxmin,zxmax,zymin,zymax,0.0d0,0.0d0)
    if (nselected .GT. 0) then 
        do nn=1,nselected
            do kk=1,surface_adtree%tree(node_select(nn))%nentry

                !Face 
                ftgt = surface_adtree%tree(node_select(nn))%entry(kk)

                !If not excluded faces 
                if ((ftgt .NE. fx1) .AND. (ftgt .NE. fx2)) then

                    !Verticies of the ends of the surface face 
                    vtx1f(:) = verticesC(surface_mesh%faces(ftgt,1),:)
                    vtx2f(:) = verticesC(surface_mesh%faces(ftgt,2),:)

                    !If an intersection then set bias and exit 
                    int_type = seg_seg_intersect_bool(vtx1,vtx2,vtx1f,vtx2f)
                    if (int_type .NE. 0.0d0) then 
                        even_bias(vv) = 1.0d0
                        exit
                    end if 
                end if 
            end do
            if (even_bias(vv) == 1.0d0) then 
                exit
            end if 
        end do 
    else
        even_bias(vv) = 0.0d0 
    end if 
end do 
do vv=1,surface_mesh%nvtx !Add from concave regions 
    if (even_bias(vv) == 0.0d0) then 

        !Next and previous vertices
        vp = get_next_vertex(surface_mesh,vv)
        vn = get_previous_vertex(surface_mesh,vv)
        
        !Next and previous vectors
        vecp(:) = verticesC(vp,:) - verticesC(vv,:)
        vecp(:) = vecp(:)/(norm2(vecp(:)) + stabval)
        vecn(:) = verticesC(vn,:) - verticesC(vv,:)
        vecn(:) = vecn(:)/(norm2(vecn(:)) + stabval)

        !dv value
        dv = dot_product(vecp,normals(vv,:)) + dot_product(vecn,normals(vv,:))

        !Add if concave 
        if (dv .GT. 0.0d0) then 
            even_bias(vv) = cm2dopt%inflayer_ew*abs(dv) + cm2dopt%inflayer_ebcbase
        end if
    end if  
end do 

!Add single encapsulated points 
do vv=1,surface_mesh%nvtx
    if (even_bias(vv) == 0.0d0) then 
        vp = get_next_vertex(surface_mesh,vv)
        vn = get_previous_vertex(surface_mesh,vv)
        if ((even_bias(vn) .NE. 0.0d0) .AND. (even_bias(vp) .NE. 0.0d0)) then 
            even_biasF(vv) = -0.5d0*(even_bias(vn) + even_bias(vp))
        else    
            even_biasF(vv) = even_bias(vv)
        end if
    else
        even_biasF(vv) = even_bias(vv)
    end if 
end do 
even_bias(:) = abs(even_biasF(:))

!Set fixed values 
ebias_fixed(:) = 0
do vv=1,surface_mesh%nvtx
    if (even_bias(vv) .NE. 0.0d0) then 
    ! if (even_bias(vv) == 1.0d0) then 
        ebias_fixed(vv) = 1
    end if
end do 

!Flood even bias 
do ff=1,cm2dopt%inflayer_enflood
    do vv=1,surface_mesh%nvtx
        if (ebias_fixed(vv) == 1) then 
            even_biasF(vv) = even_bias(vv)
        else
            vp = get_next_vertex(surface_mesh,vv)
            vn = get_previous_vertex(surface_mesh,vv)
            even_biasF(vv) = 0.5d0*(even_bias(vp) + even_bias(vn))
        end if 
    end do 
    even_bias = even_biasF
end do 

!Debug -----
! open(11,file='io/evenbvtx')
! do vv=1,surface_mesh%nvtx
!     if (even_bias(vv) .GT. 0.0d0) then 
!         write(11,*) verticesC(vv,:)
!     end if 
! end do 
! close(11)
return 
end subroutine construct_even_bias




!Evaluate growth direction subroutine ================================================
subroutine get_growth_direction(grow_dir,verticesC,even_bias,surface_mesh,steplen,cm2dopt)
implicit none 

!System data
type(cm2d_options) :: cm2dopt
type(surface_data) :: surface_mesh

!Variables - Import
real(dp) :: steplen
real(dp) :: even_bias(surface_mesh%nvtx),verticesC(surface_mesh%nvtx,2)
real(dp), dimension(:,:), allocatable :: grow_dir

!Variables - Local 
logical :: intersectstate
integer(in) :: ii,vv,nn,ff
integer(in) :: max_int_iter,vp,vn
real(dp) :: even_step,dv_step,stepstab,alpha,normvp,normvn,dnorm,lenrat,srate,h_stepav,h_total,h_totalP
real(dp) :: vecp(2),vecn(2),vtx_even(2),vec_even(2)
real(dp) :: dv(surface_mesh%nvtx),stepnorm(surface_mesh%nvtx)
real(dp) :: verticesT(surface_mesh%nvtx,2),verticesP(surface_mesh%nvtx,2),normalstep(surface_mesh%nvtx,2)
real(dp), dimension(:,:), allocatable :: normals

!Set maximum intergration iterations
max_int_iter = cm2dopt%inflayer_nintstep_max

!Set evening stepsize
even_step = cm2dopt%inflayer_ew

!Set dv step perturbation 
dv_step = cm2dopt%inflayer_dvw

!Set stabilisation parameter
stepstab = 1e-8

!Construct local normals 
call evaluate_surface_normals(normals,surface_mesh,verticesC)
if (cm2dopt%meshinout == 'in') then 
    normals(:,:) = -normals(:,:) 
end if 

!Integrate surface 
h_total = 0.0d0 
verticesP = verticesC
do ii=1,max_int_iter

    !Evaluate dv for each vertex 
    dv(:) = 0.0d0 
    do vv=1,surface_mesh%nvtx

        !Next and previous vertices
        vp = get_next_vertex(surface_mesh,vv)
        vn = get_previous_vertex(surface_mesh,vv)
        
        !Next and previous vectors
        vecp(:) = verticesP(vp,:) - verticesP(vv,:)
        vecp(:) = vecp(:)/(norm2(vecp(:)) + stepstab)
        vecn(:) = verticesP(vn,:) - verticesP(vv,:)
        vecn(:) = vecn(:)/(norm2(vecn(:)) + stepstab)

        !dv value
        dv(vv) = dot_product(vecp,normals(vv,:)) + dot_product(vecn,normals(vv,:))
    end do 
    do vv=1,surface_mesh%nvtx
        if (dv(vv) .LT. 0.0d0) then 
            dv(vv) = cm2dopt%inflayer_cvxdp*dv(vv)
        end if 
    end do 

    !Find safe step size with no self intersections 
    alpha = 1.0d0 
    do nn=1,cm2dopt%inflayer_nbclinesearch

        !Test step
        do vv=1,surface_mesh%nvtx
            srate = 1.0d0 + dv_step*dv(vv)
            verticesT(vv,1) = verticesP(vv,1) + alpha*steplen*normals(vv,1)*srate
            verticesT(vv,2) = verticesP(vv,2) + alpha*steplen*normals(vv,2)*srate
        end do 

        !Check for intersecting step vectors 
        intersectstate = has_step_intersection(verticesT,verticesP,normals,surface_mesh,cm2dopt) 

        !Shrink step if intersecting 
        if (intersectstate) then 
            alpha = 0.5d0*alpha
        else
            exit 
        end if 
    end do 
    !print *, 'alpha / nls = ',alpha,' / ',nn

    !Approximate normal step magnitude 
    do vv=1,surface_mesh%nvtx
        normalstep(vv,:) = dot_product(verticesT(vv,:) - verticesP(vv,:),normals(vv,:))*normals(vv,:)
        stepnorm(vv) = norm2(normalstep(vv,:))
    end do 

    !Find average step length 
    h_stepav = sum(stepnorm(:))/real(surface_mesh%nvtx,dp)

    !Predict total layer height 
    h_totalP = h_total + h_stepav

    !Damp step to maximum layer height if required
    alpha = 1.0d0 
    if (h_totalP .GE. steplen) then 
        alpha = abs((steplen - h_total)/h_stepav)
    end if 

    !Set new vertex positions
    verticesP(:,:) = verticesP(:,:) + alpha*(verticesT(:,:) - verticesP(:,:))

    !Accumulate approximate layer height
    h_total = h_total + alpha*h_stepav

    !Even along surface 
    do ff=1,cm2dopt%inflayer_ensubiter
        do vv=1,surface_mesh%nvtx
            if ((surface_mesh%vtx_sharp(vv) == 0) .AND. (even_bias(vv) .NE. 0.0d0)) then 

                !Next and previous vertices
                vp = get_next_vertex(surface_mesh,vv)
                vn = get_previous_vertex(surface_mesh,vv)

                !Adjacent vertex vectors
                vecp(:) = verticesP(vp,:) - verticesP(vv,:)
                vecn(:) = verticesP(vn,:) - verticesP(vv,:)

                !Evening norms
                normvp = norm2(vecp)
                normvn = norm2(vecn)
                dnorm = normvp + normvn + stepstab

                !Evening ratio
                lenrat = max((normvp/(normvn + stepstab)),(normvn/(normvp + stepstab)))

                !If large enough difference 
                if (lenrat .GE. cm2dopt%inflayer_lreb) then 

                    !Base evening vector 
                    vtx_even(:) = verticesP(vp,:)*(normvp/dnorm) + verticesP(vn,:)*(normvn/dnorm)
                    vec_even(:) = (vtx_even(:) - verticesP(vv,:)) 

                    !Help prevent shrinkage 
                    if (dot_product(vec_even(:),normals(vv,:)) .LT. 0.0d0) then 
                        vec_even(:) = vec_even(:) + cm2dopt%inflayer_enormw*norm2(vec_even(:))*normals(vv,:)
                    elseif (dot_product(vec_even(:),normals(vv,:)) .GT. 0.0d0) then 
                        vec_even(:) = vec_even(:) - cm2dopt%inflayer_enormw*norm2(vec_even(:))*normals(vv,:)
                    end if 

                    !Step
                    verticesT(vv,:) = verticesP(vv,:) + even_bias(vv)*even_step*vec_even(:) 
                else
                    verticesT(vv,:) = verticesP(vv,:)
                end if 
            else
                verticesT(vv,:) = verticesP(vv,:)
            end if 
        end do 
        verticesP = verticesT
    end do 

    !Update surface normals 
    call evaluate_surface_normals(normals,surface_mesh,verticesC)
    if (cm2dopt%meshinout == 'in') then 
        normals(:,:) = -normals(:,:) 
    end if 

    !Exit if approximate layer height has been reached by h_total
    if (h_total .GE. steplen) then 
        if (cm2dopt%dispt == 1) then
            write(*,'(A,I0,A,A,A)') '    {complete at niter = ',ii,' with height = ',real2F0_Xstring(h_total,8_in),'}'
        end if 
        exit
    end if 
end do 

!Evaluate overall layer growth direction 
grow_dir(:,:) = verticesP(:,:) - verticesC(:,:)
return 
end subroutine get_growth_direction




!Check direction intersection function ================================================
function has_step_intersection(nposN,npos0,ndir,surface_mesh,cm2dopt) result(intersectstate)
implicit none 

!Result
logical :: intersectstate

!Variables - Import
real(dp), dimension(:,:) :: nposN,npos0,ndir
type(surface_data) :: surface_mesh
type(cm2d_options) :: cm2dopt
   
!Variables - Local 
integer(in) :: vv,aa 
integer(in) :: nv,pv,cv,int_type,nacheck
real(dp) :: adjdirL
real(dp) :: lineC1(2),lineC2(2),lineA1(2),lineA2(2)

!Number of adjacent vertices to check 
nacheck = cm2dopt%inflayer_stepnacheck

!Set adjacent direction length 
adjdirL = 10.0d0 

!Set initial state
intersectstate = .false.

!Check for intersections of the steps between adjacent vertices 
do vv=1,surface_mesh%nvtx

    !Next base and previous vertices
    pv = get_previous_vertex(surface_mesh,vv)
    nv = get_next_vertex(surface_mesh,vv)

    !Current vector 
    lineC1(:) = npos0(vv,:)
    lineC2(:) = nposN(vv,:)

    !Test next
    cv = vv 
    do aa=1,nacheck

        !Test normal direction
        lineA1(:) = npos0(nv,:)
        lineA2(:) = npos0(nv,:) + adjdirL*ndir(nv,:) 
        int_type = seg_seg_intersect_bool(lineC1,lineC2,lineA1,lineA2)
        if (int_type .NE. 0) then 
            intersectstate = .true.
            return 
        end if 

        !Test segment 
        if (aa .GE. 2) then 
            lineA1(:) = npos0(cv,:)
            lineA2(:) = npos0(nv,:)
            int_type = seg_seg_intersect_bool(lineC1,lineC2,lineA1,lineA2)
            if (int_type .NE. 0) then 
                intersectstate = .true.
                return 
            end if 
        end if 

        !Update vertices 
        cv = nv
        nv = get_next_vertex(surface_mesh,nv)
    end do 

    !Test previous
    cv = vv 
    do aa=1,nacheck

        !Test normal direction
        lineA1(:) = npos0(pv,:)
        lineA2(:) = npos0(pv,:) + adjdirL*ndir(pv,:) 
        int_type = seg_seg_intersect_bool(lineC1,lineC2,lineA1,lineA2)
        if (int_type .NE. 0) then 
            intersectstate = .true.
            return 
        end if 

        !Test segment
        if (aa .GE. 2) then 
            lineA1(:) = npos0(cv,:)
            lineA2(:) = npos0(pv,:)
            int_type = seg_seg_intersect_bool(lineC1,lineC2,lineA1,lineA2)
            if (int_type .NE. 0) then 
                intersectstate = .true.
                return 
            end if 
        end if 

        !Update vertices 
        cv = pv
        pv = get_previous_vertex(surface_mesh,pv)
    end do 
end do 
return 
end function has_step_intersection




!Build inflation layer mesh subroutine ================================================
subroutine mesh_inflation_layer(volume_mesh,surface_mesh,cm2dopt)
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: vv,ilayer
integer(in) :: Nsurfedge,Nsurfvtx,Ncell_new,Nedge_new,Nvertex_new,Nlayer,lidx,fidx
integer(in) :: surf_linkindex_temp(volume_mesh%nvtx),cell_level_temp(volume_mesh%ncell)
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
cell_level_temp(:) = volume_mesh%cell_level(:)
deallocate(volume_mesh%cell_level)
allocate(volume_mesh%cell_level(Ncell_new))
volume_mesh%cell_level(1:volume_mesh%ncell) = cell_level_temp(:) 
volume_mesh%cell_level(volume_mesh%ncell:Ncell_new) = maxval(volume_mesh%cell_level(1:volume_mesh%ncell))
volume_mesh%ncell = Ncell_new

!Reset the surface mesh vertices 
surface_mesh%vertices = surface_mesh%vertices0

!Update surface link structure 
deallocate(volume_mesh%vtx_surfseg)
allocate(volume_mesh%vtx_surfseg(Nvertex_new))
volume_mesh%vtx_surfseg(:) = 0 
do vv=1,Nvertex_new

    !Link index of this vertex
    lidx = volume_mesh%surf_linkindex(vv)

    !If there is a link
    if (lidx .GT. 0) then 

        !Surface face it links to 
        fidx = volume_mesh%surf_vtx_seg(lidx)

        !Update surfseg link 
        volume_mesh%vtx_surfseg(vv) = fidx
    end if 
end do 
call build_surface_links(volume_mesh,surface_mesh)

!Clean mesh short edges
call clean_mesh_shortE(volume_mesh,cm2dopt,cm2dopt%EminLength) 
call clean_mesh_shortE(volume_mesh,cm2dopt,cm2dopt%EminLength,-1_in)

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


end module cellmesh2d_inflation_layer_mod