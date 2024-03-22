!cell_mesh2d exact mesh building module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 4.1
!Updated 22-03-2024

!Module 
module cellmesh2d_mesh_build_exact_mod
use cellmesh2d_adtree_mod
use cellmesh2d_geometry_mod
contains 

!Volume mesh construction subroutine ('exact' surface format) ===========================
subroutine construct_mesh_EXACT(volume_mesh,surface_mesh,cell_keep,vtx_external,vtx_type1,qt_mesh,surface_adtree,cm2dopt)
implicit none 

!Variables - Import
integer(in), dimension(:) :: cell_keep,vtx_external,vtx_type1
type(cm2d_options) :: cm2dopt
type(quadtree_data) :: qt_mesh
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh
type(tree_data) :: surface_adtree

!Variables - Local 
integer(in) :: type_tgt
integer(in) :: eopposite(4),celledges(4,2)
integer(in) :: cell_2_qtcell(qt_mesh%cins-1)
type(vol_mesh_data) :: volume_mesh_full

!Define opposite edges
eopposite(1) = 3
eopposite(2) = 4
eopposite(3) = 1
eopposite(4) = 2

!Define cell edges 
celledges(1,1) = 1
celledges(1,2) = 2
celledges(2,1) = 2
celledges(2,2) = 3
celledges(3,1) = 3
celledges(3,2) = 4
celledges(4,1) = 4
celledges(4,2) = 1

!Set target edge/face type for mesh 
if (cm2dopt%meshinout == 'in') then !mesh internal 
    type_tgt = 3
elseif (cm2dopt%meshinout == 'out') then !mesh external 
    type_tgt = 2
end if 

!Display
if (cm2dopt%dispt) then
    write(*,'(A)') '    {constructing base complete volume mesh}'
end if 

!Construct initial full mesh 
call build_full_mesh_EXACT(volume_mesh_full,qt_mesh,cell_keep,cell_2_qtcell,eopposite,celledges)

!Display
if (cm2dopt%dispt) then
    write(*,'(A)') '    {intersecting mesh with surface geometry}'
end if 

!Intersect the full volume mesh with the surface geometry to construct the final mesh 
call clip_vmesh2surface(volume_mesh,volume_mesh_full,surface_mesh,surface_adtree,vtx_external,vtx_type1,cm2dopt,type_tgt)
return 
end subroutine construct_mesh_EXACT




!Subroutine to clip the volume mesh to the geometry surface to construct the full volume mesh ===========================
subroutine clip_vmesh2surface(volume_mesh,volume_mesh_full,surface_mesh,surface_adtree,vtx_external,vtx_type1,cm2dopt,type_tgt)
implicit none 

!Variables - Import
integer(in) :: type_tgt
integer(in), dimension(:) :: vtx_external,vtx_type1
type(cm2d_options) :: cm2dopt
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh_full,volume_mesh
type(tree_data) :: surface_adtree

!Variables - Local
integer(in) :: ee,ii,kk,aa
integer(in) :: v1,v2,estgt,vtx_idx_smesh,int_ateend,ev1,ev2,vmedge,edir,c3,c4,intselect
integer(in) :: int_selected(cm2dopt%NintEmax),vtx_idx_odr(cm2dopt%NintEmax),vmedge_odr(cm2dopt%NintEmax)
integer(in) :: vtx_idx_vmedge(surface_mesh%nvtx)
real(dp) :: enorm
real(dp) :: vp(2),ve1(2),ve2(2),evm1(2),evm2(2),vmenormal(2),smedir(2),edge_intersect_norm(cm2dopt%NintEmax)
real(dp) :: intloc_odr(cm2dopt%NintEmax,2)
type(edgeint), dimension(:), allocatable :: edge_clip_vm,edge_clip_sm

!Build volume mesh V2E
allocate(volume_mesh_full%V2E(volume_mesh_full%nvtx,4))
volume_mesh_full%V2E(:,:) = 0 
do ee=1,volume_mesh_full%nedge
    v1 = volume_mesh_full%edge(ee,1)
    v2 = volume_mesh_full%edge(ee,2)
    do ii=1,4
        if (volume_mesh_full%V2E(v1,ii) == 0) then 
            volume_mesh_full%V2E(v1,ii) = ee 
            exit 
        end if
    end do 
    do ii=1,4
        if (volume_mesh_full%V2E(v2,ii) == 0) then 
            volume_mesh_full%V2E(v2,ii) = ee 
            exit 
        end if
    end do 
end do

!Initialise clipped surface structures 
vtx_idx_vmedge(:) = 0 
allocate(edge_clip_sm(surface_mesh%nfcs))
do ee=1,surface_mesh%nfcs
    edge_clip_sm(ee)%nint = 0 
    edge_clip_sm(ee)%type = 0 
    allocate(edge_clip_sm(ee)%intloc(cm2dopt%NintEmax,2))
    allocate(edge_clip_sm(ee)%vtx_idx(cm2dopt%NintEmax))
    allocate(edge_clip_sm(ee)%vmedge(cm2dopt%NintEmax))
    allocate(edge_clip_sm(ee)%vmeshcell(cm2dopt%NintEmax))
    edge_clip_sm(ee)%vtx_idx(:) = 0 
    edge_clip_sm(ee)%vmedge(:) = 0 
    edge_clip_sm(ee)%vmeshcell(:) = 0 
    edge_clip_sm(ee)%refend1 = surface_mesh%faces(ee,1)
end do 

!Set index for surface mesh intersection vertices 
vtx_idx_smesh = surface_mesh%nvtx

!Intersect volume mesh edges with the surface mesh
call clip_edges2surface(edge_clip_vm,vtx_external,vtx_type1,volume_mesh_full,surface_mesh,surface_adtree,cm2dopt)

!Accumulate intersect vertices to the clipped edge surface mesh structure 
do ee=1,volume_mesh_full%nedge
    if (edge_clip_vm(ee)%nint .NE. 0) then 
        do ii=1,edge_clip_vm(ee)%nint

            !Target surface edge 
            estgt = edge_clip_vm(ee)%surfseg(ii)

            !Edge ends
            ve1(:) = surface_mesh%vertices(surface_mesh%faces(estgt,1),:)
            ve2(:) = surface_mesh%vertices(surface_mesh%faces(estgt,2),:)

            !Intersection location 
            vp(:) = edge_clip_vm(ee)%intloc(ii,:)

            !Edge norm 
            enorm = norm2(ve2(:) - ve1(:))

            !Check if within tollerance of an end of this edge
            int_ateend = 0 
            if (norm2(ve2(:) - vp(:)) .LE. set_cmp_prc(enorm*cm2dopt%intcointol,min_precision)) then 
                int_ateend = surface_mesh%faces(estgt,2)
                vp(:) = surface_mesh%vertices(surface_mesh%faces(estgt,2),:)
            elseif (norm2(ve1(:) - vp(:)) .LE. set_cmp_prc(enorm*cm2dopt%intcointol,min_precision)) then 
                int_ateend = surface_mesh%faces(estgt,1)
                vp(:) = surface_mesh%vertices(surface_mesh%faces(estgt,1),:)
            end if

            !Set intersection correctly as edge or vertex intersect on the surface mesh 
            if (int_ateend .NE. 0) then !Set as vertex intersect

                !Set intersect index on the volume mesh edge as the surface vertex it intersects
                edge_clip_vm(ee)%vtx_idx(ii) = int_ateend !set as positive as this points to the same vertex set as the surface mesh vertices ?? ***********
                    
                !Tag this surface vertex with this volume mesh edge 
                vtx_idx_vmedge(int_ateend) = ee 
            else !Set as edge intersect 

                !Increment vertex index
                vtx_idx_smesh = vtx_idx_smesh + 1 

                !Set properties 
                edge_clip_sm(estgt)%type = 4
                edge_clip_sm(estgt)%nint = edge_clip_sm(estgt)%nint + 1
                edge_clip_sm(estgt)%intloc(edge_clip_sm(estgt)%nint,:) = edge_clip_vm(ee)%intloc(ii,:)
                edge_clip_sm(estgt)%vmedge(edge_clip_sm(estgt)%nint) = ee 

                !Store this vertex index
                edge_clip_sm(estgt)%vtx_idx(edge_clip_sm(estgt)%nint) = vtx_idx_smesh !edge_clip_vm(ee)%vtx_idx(ii)
                edge_clip_vm(ee)%vtx_idx(ii) = vtx_idx_smesh
            end if 
        end do 
    end if
end do 

!Order surface mesh edge intersections from the edge start 
do ee=1,surface_mesh%nfcs
    if (edge_clip_sm(ee)%type == 4) then 

        !Edge ends 
        ve1(:) = surface_mesh%vertices(surface_mesh%faces(ee,1),:)
        ve2(:) = surface_mesh%vertices(surface_mesh%faces(ee,2),:)

        !Distances
        edge_intersect_norm(:) = 0.0d0  
        do kk=1,edge_clip_sm(ee)%nint
            edge_intersect_norm(kk) = norm2(edge_clip_sm(ee)%intloc(kk,:) - ve1(:))
        end do 

        !Order
        int_selected(:) = 0 
        do kk=1,edge_clip_sm(ee)%nint

            !Select intersection at maximum distance
            intselect = maxloc(edge_intersect_norm(1:edge_clip_sm(ee)%nint),1,&
                               mask = int_selected(1:edge_clip_sm(ee)%nint) == 0)

            !Store 
            intloc_odr(edge_clip_sm(ee)%nint-kk+1,:) = edge_clip_sm(ee)%intloc(intselect,:)
            vtx_idx_odr(edge_clip_sm(ee)%nint-kk+1) = edge_clip_sm(ee)%vtx_idx(intselect)
            vmedge_odr(edge_clip_sm(ee)%nint-kk+1) = edge_clip_sm(ee)%vmedge(intselect)

            !Set distance to zero and set to visited
            int_selected(intselect) = 1
            edge_intersect_norm(intselect) = 0.0d0 
        end do 

        !Store
        edge_clip_sm(ee)%intloc(1:edge_clip_sm(ee)%nint,:) = intloc_odr(1:edge_clip_sm(ee)%nint,:)
        edge_clip_sm(ee)%vtx_idx(1:edge_clip_sm(ee)%nint) = vtx_idx_odr(1:edge_clip_sm(ee)%nint)
        edge_clip_sm(ee)%vmedge(1:edge_clip_sm(ee)%nint) = vmedge_odr(1:edge_clip_sm(ee)%nint)
    end if 
end do 

!Build clipped sub-edge mesh on clipped surface mesh edges and set the volume mesh cell for each segment 
do ee=1,surface_mesh%nfcs

    !Surface mesh edge ends 
    ev1 = surface_mesh%faces(ee,1)
    ev2 = surface_mesh%faces(ee,2)

    !Surface mesh edge direction 
    smedir(:) = surface_mesh%vertices(ev2,:) - surface_mesh%vertices(ev1,:)
    
    !Construct
    if (edge_clip_sm(ee)%type == 4) then 

        !Allocate
        allocate(edge_clip_sm(ee)%edge_mesh(edge_clip_sm(ee)%nint+1,2))
                    
        !Set first segment (negative numbers index the intersections along the edge)
        edge_clip_sm(ee)%edge_mesh(1,1) = ev1
        edge_clip_sm(ee)%edge_mesh(1,2) = -1

        !Set middle segments 
        if (edge_clip_sm(ee)%nint .GE. 2) then 
            do kk=2,edge_clip_sm(ee)%nint
                edge_clip_sm(ee)%edge_mesh(kk,1) = -1*(kk - 1)
                edge_clip_sm(ee)%edge_mesh(kk,2) = -1*kk
            end do 
        end if 

        !Set final segment 
        edge_clip_sm(ee)%edge_mesh(edge_clip_sm(ee)%nint+1,1) = -edge_clip_sm(ee)%nint
        edge_clip_sm(ee)%edge_mesh(edge_clip_sm(ee)%nint+1,2) = ev2

        !For each sub-edge segment assign a volume mesh cell 
        do ii=1,edge_clip_sm(ee)%nint !Set from edge intersections

            !Volume mesh edge on this intersect
            vmedge = edge_clip_sm(ee)%vmedge(ii)

            !Cells on this edge 
            c3 = volume_mesh_full%edge(vmedge,3)
            c4 = volume_mesh_full%edge(vmedge,4)

            !Volume mesh edge direction and normal 
            evm1(:) = volume_mesh_full%vertices(volume_mesh_full%edge(vmedge,1),:)
            evm2(:) = volume_mesh_full%vertices(volume_mesh_full%edge(vmedge,2),:)
            vmenormal(1) = evm1(2) - evm2(2)
            vmenormal(2) = evm2(1) - evm1(1)

            !Set direction 
            if (dot_product(smedir,vmenormal) .GT. 0.0d0) then 
                edir = 1
            elseif (dot_product(smedir,vmenormal) .LT. 0.0d0) then 
                edir = -1
            else
                edir = 0
                print *, '** indeterminate surface edge - volume edge cross direction at edges (s/v) ',ee,vmedge
                ! print *, evm1(:)
                ! print *, evm2(:)
            end if

            !Add cells to all sub-edges on this vertex with intersections
            do aa=1,edge_clip_sm(ee)%nint+1
                if (edge_clip_sm(ee)%vmeshcell(aa) == 0) then !this sub-edge has no assigned volume mesh cell
                    v1 = edge_clip_sm(ee)%edge_mesh(aa,1)
                    v2 = edge_clip_sm(ee)%edge_mesh(aa,2)
                    if (v1 .LT. 0) then 
                        if (abs(v1) == ii) then !sub-edge starts at this intersection
                            if (edir == 1) then 
                                edge_clip_sm(ee)%vmeshcell(aa) = c3 
                            elseif (edir == -1) then 
                                edge_clip_sm(ee)%vmeshcell(aa) = c4 
                            end if 
                        end if 
                    end if 
                    if (v2 .LT. 0) then 
                        if (abs(v2) == ii) then !sub-edge ends at this intersection
                            if (edir == 1) then 
                                edge_clip_sm(ee)%vmeshcell(aa) = c4 
                            elseif (edir == -1) then 
                                edge_clip_sm(ee)%vmeshcell(aa) = c3 
                            end if 
                        end if 
                    end if 
                end if
            end do 
        end do 

        !Set from vertex intersections 
        if (vtx_idx_vmedge(ev1) .NE. 0) then !end 1 (first segment)

            !Volume mesh edge on this intersect
            vmedge = vtx_idx_vmedge(ev1)

            !Cells on this edge 
            c3 = volume_mesh_full%edge(vmedge,3)
            c4 = volume_mesh_full%edge(vmedge,4)

            !Volume mesh edge direction and normal 
            evm1(:) = volume_mesh_full%vertices(volume_mesh_full%edge(vmedge,1),:)
            evm2(:) = volume_mesh_full%vertices(volume_mesh_full%edge(vmedge,2),:)
            vmenormal(1) = evm1(2) - evm2(2)
            vmenormal(2) = evm2(1) - evm1(1)

            !Set direction 
            if (dot_product(smedir,vmenormal) .GT. 0.0d0) then 
                edir = 1
            elseif (dot_product(smedir,vmenormal) .LT. 0.0d0) then 
                edir = -1
            else
                edir = 0
                print *, '** indeterminate surface edge - volume edge cross direction at edges (s/v) ',ee,vmedge
                ! print *, evm1(:)
                ! print *, evm2(:)
            end if

            !Set cell of first segment 
            if (edge_clip_sm(ee)%vmeshcell(1) == 0) then
                if (edir == 1) then 
                    edge_clip_sm(ee)%vmeshcell(1) = c3 
                elseif (edir == -1) then 
                    edge_clip_sm(ee)%vmeshcell(1) = c4 
                end if 
            end if 
        end if
        if (vtx_idx_vmedge(ev2) .NE. 0) then !end 2 (final segment)

            !Volume mesh edge on this intersect
            vmedge = vtx_idx_vmedge(ev2)

            !Cells on this edge 
            c3 = volume_mesh_full%edge(vmedge,3)
            c4 = volume_mesh_full%edge(vmedge,4)

            !Volume mesh edge direction and normal 
            evm1(:) = volume_mesh_full%vertices(volume_mesh_full%edge(vmedge,1),:)
            evm2(:) = volume_mesh_full%vertices(volume_mesh_full%edge(vmedge,2),:)
            vmenormal(1) = evm1(2) - evm2(2)
            vmenormal(2) = evm2(1) - evm1(1)

            !Set direction 
            if (dot_product(smedir,vmenormal) .GT. 0.0d0) then 
                edir = 1
            elseif (dot_product(smedir,vmenormal) .LT. 0.0d0) then 
                edir = -1
            else
                edir = 0
                print *, '** indeterminate surface edge - volume edge cross direction at edges (s/v) ',ee,vmedge
                ! print *, evm1(:)
                ! print *, evm2(:)
            end if

            !Set cell of final segment 
            if (edge_clip_sm(ee)%vmeshcell(edge_clip_sm(ee)%nint+1) == 0) then
                if (edir == 1) then 
                    edge_clip_sm(ee)%vmeshcell(edge_clip_sm(ee)%nint+1) = c4 
                elseif (edir == -1) then 
                    edge_clip_sm(ee)%vmeshcell(edge_clip_sm(ee)%nint+1) = c3 
                end if 
            end if 
        end if 
    else !Non intersecting edge -> set volume mesh cell if either end is a vertex intersect
        if (vtx_idx_vmedge(ev1) .NE. 0) then !end 1 (first segment)

            !Volume mesh edge on this intersect
            vmedge = vtx_idx_vmedge(ev1)

            !Cells on this edge 
            c3 = volume_mesh_full%edge(vmedge,3)
            c4 = volume_mesh_full%edge(vmedge,4)

            !Volume mesh edge direction and normal 
            evm1(:) = volume_mesh_full%vertices(volume_mesh_full%edge(vmedge,1),:)
            evm2(:) = volume_mesh_full%vertices(volume_mesh_full%edge(vmedge,2),:)
            vmenormal(1) = evm1(2) - evm2(2)
            vmenormal(2) = evm2(1) - evm1(1)

            !Set direction 
            if (dot_product(smedir,vmenormal) .GT. 0.0d0) then 
                edir = 1
            elseif (dot_product(smedir,vmenormal) .LT. 0.0d0) then 
                edir = -1
            else
                edir = 0
                print *, '** indeterminate surface edge - volume edge cross direction at edges (s/v) ',ee,vmedge
                ! print *, evm1(:)
                ! print *, evm2(:)
            end if

            !Set cell of segment 
            if (edge_clip_sm(ee)%vmeshcell(1) == 0) then
                if (edir == 1) then 
                    edge_clip_sm(ee)%vmeshcell(1) = c3 
                elseif (edir == -1) then 
                    edge_clip_sm(ee)%vmeshcell(1) = c4 
                end if 
            end if  
        end if
        if (vtx_idx_vmedge(ev2) .NE. 0) then !end 2

            !Volume mesh edge on this intersect
            vmedge = vtx_idx_vmedge(ev2)

            !Cells on this edge 
            c3 = volume_mesh_full%edge(vmedge,3)
            c4 = volume_mesh_full%edge(vmedge,4)

            !Volume mesh edge direction and normal 
            evm1(:) = volume_mesh_full%vertices(volume_mesh_full%edge(vmedge,1),:)
            evm2(:) = volume_mesh_full%vertices(volume_mesh_full%edge(vmedge,2),:)
            vmenormal(1) = evm1(2) - evm2(2)
            vmenormal(2) = evm2(1) - evm1(1)

            !Set direction 
            if (dot_product(smedir,vmenormal) .GT. 0.0d0) then 
                edir = 1
            elseif (dot_product(smedir,vmenormal) .LT. 0.0d0) then 
                edir = -1
            else
                edir = 0
                print *, '** indeterminate surface edge - volume edge cross direction at edges (s/v) ',ee,vmedge
                print *, evm1(:)
                print *, evm2(:)
            end if

            !Set cell of segment 
            if (edge_clip_sm(ee)%vmeshcell(1) == 0) then
                if (edir == 1) then 
                    edge_clip_sm(ee)%vmeshcell(1) = c4 
                elseif (edir == -1) then 
                    edge_clip_sm(ee)%vmeshcell(1) = c3 
                end if 
            end if  
        end if 
    end if 
end do 

!Build complete list of surface vertices including all volume mesh to surface mesh clipping intersections 
surface_mesh%nvtxf = vtx_idx_smesh
allocate(surface_mesh%vertices_full(vtx_idx_smesh,2))
surface_mesh%vertices_full(1:surface_mesh%nvtx,:) = surface_mesh%vertices(:,:)
do ee=1,surface_mesh%nfcs
    if (edge_clip_sm(ee)%nint .GT. 0) then 
        do ii=1,edge_clip_sm(ee)%nint
            surface_mesh%vertices_full(edge_clip_sm(ee)%vtx_idx(ii),:) = edge_clip_sm(ee)%intloc(ii,:)
        end do
    end if
end do 

!Construct complete final volume mesh 
call build_final_vmesh(volume_mesh,volume_mesh_full,surface_mesh,edge_clip_vm,edge_clip_sm,vtx_idx_vmedge,&
                       cm2dopt,type_tgt)
return
end subroutine clip_vmesh2surface




!Build final clipped volume mesh subroutine ===========================
subroutine build_final_vmesh(volume_mesh,volume_mesh_full,surface_mesh,edge_clip_vm,edge_clip_sm,vtx_idx_vmedge,&
                             cm2dopt,type_tgt)
implicit none 

!System data
type(edgeint), dimension(:) :: edge_clip_vm,edge_clip_sm
type(vol_mesh_data) :: volume_mesh_full,volume_mesh
type(surface_data) :: surface_mesh
type(cm2d_options) :: cm2dopt

!Variables - Import
integer(in) :: type_tgt
integer(in) :: vtx_idx_vmedge(surface_mesh%nvtx)

!Variables - Local
integer(in) :: ii,ee,vv,ff
integer(in) :: Nvtx_vmf,Nedge_vmf,ev1,ev2,Nesurf,v1,v2,Nvtx_surf,e1,e2,c1,c2,nupdate,l1,l2
integer(in) :: vtxmap_vm(volume_mesh_full%nvtx),vtxmap_sm(surface_mesh%nvtxf)
integer(in) :: smvtx2vmedge(surface_mesh%nvtx,2)
integer(in), dimension(:), allocatable :: vtx_surfint,surf_vertices,surf_vtx_idx
integer(in), dimension(:,:), allocatable :: edge_surf,v2e_vsurf
real(dp) :: r1,r2

!Initialise
Nedge_vmf = 0 
Nvtx_vmf = 0
vtxmap_vm(:) = 0 
vtxmap_sm(:) = 0 

!Index retained vertices in the volume mesh 
do ee=1,volume_mesh_full%nedge
    if (edge_clip_vm(ee)%type == type_tgt) then 
        ev1 = volume_mesh_full%edge(ee,1)
        ev2 = volume_mesh_full%edge(ee,2)
        if (vtxmap_vm(ev1) == 0) then 
            Nvtx_vmf = Nvtx_vmf + 1
            vtxmap_vm(ev1) = Nvtx_vmf
        end if 
        if (vtxmap_vm(ev2) == 0) then 
            Nvtx_vmf = Nvtx_vmf + 1
            vtxmap_vm(ev2) = Nvtx_vmf
        end if 
    elseif (edge_clip_vm(ee)%type == 4) then 
        do ii=1,edge_clip_vm(ee)%nint+1
            if (edge_clip_vm(ee)%int_seg_type(ii) == type_tgt) then 
                ev1 = edge_clip_vm(ee)%edge_mesh(ii,1)
                ev2 = edge_clip_vm(ee)%edge_mesh(ii,2)
                if (ev1 .GT. 0) then 
                    if (vtxmap_vm(ev1) == 0) then 
                        Nvtx_vmf = Nvtx_vmf + 1
                        vtxmap_vm(ev1) = Nvtx_vmf
                    end if 
                end if 
                if (ev2 .GT. 0) then 
                    if (vtxmap_vm(ev2) == 0) then 
                        Nvtx_vmf = Nvtx_vmf + 1
                        vtxmap_vm(ev2) = Nvtx_vmf
                    end if 
                end if 
            end if 
        end do
    end if 
end do 
! open(11,file='vtxret.dat')
! do vv=1,volume_mesh_full%nvtx
!     if (vtxmap_vm(vv) .NE. 0) then 
!         write(11,*) volume_mesh_full%vertices(vv,:)
!     end if 
! end do 
! close(11)

!Index vertices in the surface mesh 
do vv=1,surface_mesh%nvtxf
    Nvtx_vmf = Nvtx_vmf + 1
    vtxmap_sm(vv) = Nvtx_vmf 
end do 

!Build complete vertex list
volume_mesh%nvtx = Nvtx_vmf
allocate(volume_mesh%vertices(Nvtx_vmf,2))
do vv=1,volume_mesh_full%nvtx
    if (vtxmap_vm(vv) .NE. 0) then 
        volume_mesh%vertices(vtxmap_vm(vv),:) = volume_mesh_full%vertices(vv,:)
    end if 
end do 
do vv=1,surface_mesh%nvtxf
    volume_mesh%vertices(vtxmap_sm(vv),:) = surface_mesh%vertices_full(vv,:)
end do 

!Index complete edge list in volume mesh 
do ee=1,volume_mesh_full%nedge
    if (edge_clip_vm(ee)%type == type_tgt) then 
        Nedge_vmf = Nedge_vmf + 1
    elseif (edge_clip_vm(ee)%type == 4) then 
        do ii=1,edge_clip_vm(ee)%nint+1
            if (edge_clip_vm(ee)%int_seg_type(ii) == type_tgt) then 
                Nedge_vmf = Nedge_vmf + 1
            end if
        end do 
    end if 
end do 

!Index complete edge list in surface mesh 
do ee=1,surface_mesh%nfcs
    if (edge_clip_sm(ee)%type == 4) then 
        Nedge_vmf = Nedge_vmf + edge_clip_sm(ee)%nint + 1
    else
        Nedge_vmf = Nedge_vmf + 1
    end if 
end do 

!Build complete list of mesh edges 
volume_mesh%nedge = Nedge_vmf
allocate(volume_mesh%edge(Nedge_vmf,4))
volume_mesh%edge(:,:) = 0 
Nedge_vmf = 0 
do ee=1,volume_mesh_full%nedge !From volume mesh
    if (edge_clip_vm(ee)%type == type_tgt) then 
        Nedge_vmf = Nedge_vmf + 1
        volume_mesh%edge(Nedge_vmf,1) = vtxmap_vm(volume_mesh_full%edge(ee,1))
        volume_mesh%edge(Nedge_vmf,2) = vtxmap_vm(volume_mesh_full%edge(ee,2))
        volume_mesh%edge(Nedge_vmf,3) = volume_mesh_full%edge(ee,3)
        volume_mesh%edge(Nedge_vmf,4) = volume_mesh_full%edge(ee,4)
    elseif (edge_clip_vm(ee)%type == 4) then 
        do ii=1,edge_clip_vm(ee)%nint+1
            if (edge_clip_vm(ee)%int_seg_type(ii) == type_tgt) then 
                ev1 = edge_clip_vm(ee)%edge_mesh(ii,1)
                ev2 = edge_clip_vm(ee)%edge_mesh(ii,2)
                if (ev1 .GT. 0) then 
                    ev1 = vtxmap_vm(volume_mesh_full%edge(ee,1))
                else
                    ev1 = vtxmap_sm(edge_clip_vm(ee)%vtx_idx(abs(ev1)))
                end if 
                if (ev2 .GT. 0) then 
                    ev2 = vtxmap_vm(volume_mesh_full%edge(ee,2))
                else
                    ev2 = vtxmap_sm(edge_clip_vm(ee)%vtx_idx(abs(ev2)))
                end if 
                Nedge_vmf = Nedge_vmf + 1
                volume_mesh%edge(Nedge_vmf,1) = ev1
                volume_mesh%edge(Nedge_vmf,2) = ev2
                volume_mesh%edge(Nedge_vmf,3) = volume_mesh_full%edge(ee,3)
                volume_mesh%edge(Nedge_vmf,4) = volume_mesh_full%edge(ee,4)
            end if
        end do 
    end if
end do 
smvtx2vmedge(:,:) = 0 
do ee=1,surface_mesh%nfcs !From surface mesh 
    if (edge_clip_sm(ee)%type == 4) then 
        do ii=1,edge_clip_sm(ee)%nint+1
            ev1 = edge_clip_sm(ee)%edge_mesh(ii,1)
            ev2 = edge_clip_sm(ee)%edge_mesh(ii,2)
            if (ev1 .GT. 0) then 
                ev1 = vtxmap_sm(surface_mesh%faces(ee,1))
            else
                ev1 = vtxmap_sm(edge_clip_sm(ee)%vtx_idx(abs(ev1)))
            end if 
            if (ev2 .GT. 0) then 
                ev2 = vtxmap_sm(surface_mesh%faces(ee,2))
            else
                ev2 = vtxmap_sm(edge_clip_sm(ee)%vtx_idx(abs(ev2)))
            end if 
            Nedge_vmf = Nedge_vmf + 1
            if (cm2dopt%meshinout == 'in') then !mesh internal 
                volume_mesh%edge(Nedge_vmf,1) = ev1
                volume_mesh%edge(Nedge_vmf,2) = ev2
            elseif (cm2dopt%meshinout == 'out') then !mesh external 
                volume_mesh%edge(Nedge_vmf,1) = ev2
                volume_mesh%edge(Nedge_vmf,2) = ev1
            end if 
            volume_mesh%edge(Nedge_vmf,3) = -1
            volume_mesh%edge(Nedge_vmf,4) = edge_clip_sm(ee)%vmeshcell(ii)
            if (edge_clip_sm(ee)%edge_mesh(ii,1) .GT. 0) then 
                if (smvtx2vmedge(surface_mesh%faces(ee,1),1) == 0) then 
                    smvtx2vmedge(surface_mesh%faces(ee,1),1) = Nedge_vmf
                elseif (smvtx2vmedge(surface_mesh%faces(ee,1),2) == 0) then 
                    smvtx2vmedge(surface_mesh%faces(ee,1),2) = Nedge_vmf
                else
                    print *, '** surface vertex edge link overflow: ',smvtx2vmedge(surface_mesh%faces(ee,1),:),Nedge_vmf
                end if
            end if 
            if (edge_clip_sm(ee)%edge_mesh(ii,2) .GT. 0) then 
                if (smvtx2vmedge(surface_mesh%faces(ee,2),1) == 0) then 
                    smvtx2vmedge(surface_mesh%faces(ee,2),1) = Nedge_vmf
                elseif (smvtx2vmedge(surface_mesh%faces(ee,2),2) == 0) then 
                    smvtx2vmedge(surface_mesh%faces(ee,2),2) = Nedge_vmf
                else
                    print *, '** surface vertex edge link overflow: ',smvtx2vmedge(surface_mesh%faces(ee,2),:),Nedge_vmf
                end if
            end if 
        end do 
    else
        Nedge_vmf = Nedge_vmf + 1
        if (cm2dopt%meshinout == 'in') then !mesh internal 
            volume_mesh%edge(Nedge_vmf,1) = vtxmap_sm(surface_mesh%faces(ee,1))
            volume_mesh%edge(Nedge_vmf,2) = vtxmap_sm(surface_mesh%faces(ee,2))
        elseif (cm2dopt%meshinout == 'out') then !mesh external 
            volume_mesh%edge(Nedge_vmf,1) = vtxmap_sm(surface_mesh%faces(ee,2))
            volume_mesh%edge(Nedge_vmf,2) = vtxmap_sm(surface_mesh%faces(ee,1))
        end if 
        volume_mesh%edge(Nedge_vmf,3) = -1
        volume_mesh%edge(Nedge_vmf,4) = edge_clip_sm(ee)%vmeshcell(1)
        if (smvtx2vmedge(surface_mesh%faces(ee,1),1) == 0) then 
            smvtx2vmedge(surface_mesh%faces(ee,1),1) = Nedge_vmf
        elseif (smvtx2vmedge(surface_mesh%faces(ee,1),2) == 0) then 
            smvtx2vmedge(surface_mesh%faces(ee,1),2) = Nedge_vmf
        else
            print *, '** surface vertex edge link overflow: ',smvtx2vmedge(surface_mesh%faces(ee,1),:),Nedge_vmf
        end if
        if (smvtx2vmedge(surface_mesh%faces(ee,2),1) == 0) then 
            smvtx2vmedge(surface_mesh%faces(ee,2),1) = Nedge_vmf
        elseif (smvtx2vmedge(surface_mesh%faces(ee,2),2) == 0) then 
            smvtx2vmedge(surface_mesh%faces(ee,2),2) = Nedge_vmf
        else
            print *, '** surface vertex edge link overflow: ',smvtx2vmedge(surface_mesh%faces(ee,2),:),Nedge_vmf
        end if
    end if 
end do 
allocate(surface_mesh%smvtx2vmedge(surface_mesh%nvtx,2))
surface_mesh%smvtx2vmedge(:,:) = smvtx2vmedge(:,:)

!Extract surface mesh within volume mesh 
Nesurf = 0 
allocate(edge_surf(volume_mesh%nedge,3))
edge_surf(:,:) = 0
do ee=1,volume_mesh%nedge
    if ((volume_mesh%edge(ee,3) == -1) .OR. (volume_mesh%edge(ee,4) == -1)) then 
        Nesurf = Nesurf + 1
        edge_surf(Nesurf,1:2) = volume_mesh%edge(ee,1:2)
        edge_surf(Nesurf,3) = ee 
    end if 
end do 

!Build v2e for this surface sub-mesh
allocate(v2e_vsurf(volume_mesh%nvtx,2))
v2e_vsurf(:,:) = 0 
do ee=1,Nesurf
    v1 = edge_surf(ee,1)
    v2 = edge_surf(ee,2)
    v2e_vsurf(v1,2) = edge_surf(ee,3) 
    v2e_vsurf(v2,1) = edge_surf(ee,3) 
end do 

!Tag vertices that are surface intersections 
allocate(vtx_surfint(volume_mesh%nvtx))
vtx_surfint(:) = 0 
do ee=1,volume_mesh_full%nedge
    if (edge_clip_vm(ee)%type == 4) then 
        do ii=1,edge_clip_vm(ee)%nint+1
            if (edge_clip_vm(ee)%int_seg_type(ii) == type_tgt) then 
                ev1 = edge_clip_vm(ee)%edge_mesh(ii,1)
                ev2 = edge_clip_vm(ee)%edge_mesh(ii,2)
                if (ev1 .LT. 0) then 
                    ev1 = vtxmap_sm(edge_clip_vm(ee)%vtx_idx(abs(ev1)))
                    vtx_surfint(ev1) = 1
                end if 
                if (ev2 .LT. 0) then 
                    ev2 = vtxmap_sm(edge_clip_vm(ee)%vtx_idx(abs(ev2)))
                    vtx_surfint(ev2) = 1
                end if 
            end if 
        end do
    end if 
end do 
do vv=1,surface_mesh%nvtx
    if (vtx_idx_vmedge(vv) .NE. 0) then 
        vtx_surfint(vtxmap_sm(vv)) = 1
    end if 
end do 
! open(11,file='io/surfvtx')
!     do vv=1,volume_mesh%nvtx
!         if (vtx_surfint(vv) == 1) then 
!             write(11,*) volume_mesh%vertices(vv,:)
!         end if 
!     end do 
! close(11)

!Index surface vertices that are not volume mesh - surface intersections 
Nvtx_surf = 0 
allocate(surf_vtx_idx(volume_mesh%nvtx))
surf_vtx_idx(:) = 0 
do ee=1,volume_mesh%nedge
    if (volume_mesh%edge(ee,3) == -1) then 
        ev1 = volume_mesh%edge(ee,1)
        ev2 = volume_mesh%edge(ee,2)
        if (vtx_surfint(ev1) == 0) then 
            if (surf_vtx_idx(ev1) == 0) then 
                Nvtx_surf = Nvtx_surf + 1
                surf_vtx_idx(ev1) = Nvtx_surf
            end if 
        end if 
        if (vtx_surfint(ev2) == 0) then 
            if (surf_vtx_idx(ev2) == 0) then 
                Nvtx_surf = Nvtx_surf + 1
                surf_vtx_idx(ev2) = Nvtx_surf
            end if 
        end if 
    end if 
end do 
allocate(surf_vertices(Nvtx_surf))
surf_vertices(:) = 0 
do vv=1,volume_mesh%nvtx
    if (surf_vtx_idx(vv) .NE. 0) then 
        surf_vertices(surf_vtx_idx(vv)) = vv 
    end if
end do 
! open(11,file='io/surfvtx')
!     do vv=1,Nvtx_surf
!         write(11,*) volume_mesh%vertices(surf_vertices(vv),:)
!     end do 
! close(11)

!Flood cell associations of surface mesh edges along the surface 
do ff=1,Nvtx_surf
    nupdate = 0 
    do vv=1,Nvtx_surf
        e1 = v2e_vsurf(surf_vertices(vv),1)
        e2 = v2e_vsurf(surf_vertices(vv),2)
        if ((volume_mesh%edge(e1,4) == 0) .AND. (volume_mesh%edge(e2,4) .NE. 0)) then 
            volume_mesh%edge(e1,4) = volume_mesh%edge(e2,4)
            nupdate = nupdate + 1
        end if 
        if ((volume_mesh%edge(e2,4) == 0) .AND. (volume_mesh%edge(e1,4) .NE. 0)) then 
            volume_mesh%edge(e2,4) = volume_mesh%edge(e1,4)
            nupdate = nupdate + 1
        end if 
    end do 
    !print *,nupdate
    if (nupdate == 0) then !Exit at complete flood 
        exit
    end if  
end do 

!Set the initial number of cells on the volume mesh to that of the full mesh (trimmed later to the true value when cells are re-indexed)
volume_mesh%ncell = volume_mesh_full%ncell

!Set cell levels (trimmed later to the true cell indexes only when cells are re-indexed)
allocate(volume_mesh%cell_level(volume_mesh%ncell))
volume_mesh%cell_level(:) = volume_mesh_full%cell_level(:)

!Allocate surface segment links 
allocate(volume_mesh%vtx_surfseg(volume_mesh%nvtx))
volume_mesh%vtx_surfseg(:) = 0 

!Build surface segment links from mesh-surface intersections
do ee=1,volume_mesh_full%nedge
    if (edge_clip_vm(ee)%type == 4) then 
        do ii=1,edge_clip_vm(ee)%nint
            v1 = vtxmap_sm(edge_clip_vm(ee)%vtx_idx(ii))
            volume_mesh%vtx_surfseg(v1) = edge_clip_vm(ee)%surfseg(ii)
        end do
    end if 
end do 

!Build surface segment links from the base surface mesh ***


!Set cell size based search radius for each surface vertex
allocate(surface_mesh%vtx_rsearch(surface_mesh%nvtx))
surface_mesh%vtx_rsearch(:) = 0.0d0  
do vv=1,surface_mesh%nvtx

    !Edges on this vertex
    e1 = surface_mesh%smvtx2vmedge(vv,1)  
    e2 = surface_mesh%smvtx2vmedge(vv,2)  

    !Cells and levels 
    if (e1 .GT. 0) then 
        c1 = volume_mesh%edge(e1,4) 
        if (c1 .GT. 0) then 
            l1 = volume_mesh%cell_level(c1)
            r1 = 2.0d0*cm2dopt%far_field_bound/(2.0d0**(l1 - 1))
        else
            r1 = 0.0d0
        end if 
    else
        r1 = 0.0d0 
    end if
    if (e2 .GT. 0) then
        c2 = volume_mesh%edge(e1,4) 
        if (c2 .GT. 0) then 
            l2 = volume_mesh%cell_level(c2)
            r2 = 2.0d0*cm2dopt%far_field_bound/(2.0d0**(l2 - 1))
        else
            r2 = 0.0d0 
        end if 
    else
        r2 = 0.0d0 
    end if
    
    !Select largest bound to be the search radius 
    surface_mesh%vtx_rsearch(vv) = max(r1,r2)    
    !print *, surface_mesh%vtx_rsearch(vv)
end do 

!Debug -----
! print *, Nvtx_vmf,Nedge_vmf,volume_mesh%ncell
return 
end subroutine build_final_vmesh




!Clip edges to surface subroutine ===========================
subroutine clip_edges2surface(edge_clip,vtx_external,vtx_type1,volume_mesh_full,surface_mesh,surface_adtree,cm2dopt)
implicit none 

!Variables - Import
integer(in), dimension(:) :: vtx_external,vtx_type1
type(vol_mesh_data) :: volume_mesh_full
type(edgeint), dimension(:), allocatable :: edge_clip
type(cm2d_options) :: cm2dopt
type(surface_data) :: surface_mesh
type(tree_data) :: surface_adtree

!Variables - Local
integer(in) :: ee,nn,kk,ii,aa
integer(in) :: ev1,ev2,echeck,nselected,N_edge_intersect,segval,segnew,intnew,int_unique,vref,eadj,etype_check,Neintersect
integer(in) :: intselect,sint_idx
integer(in) :: node_select(surface_adtree%nnode),edge_check(volume_mesh_full%nedge),int_selected(cm2dopt%NintEmax)
integer(in) :: edge_intersect_seg(cm2dopt%NintEmax),edge_intersect_type(cm2dopt%NintEmax)
real(dp) :: segnorm,cpadSZ,zxmin,zxmax,zymin,zymax,zzmin,zzmax,efrac,idist
real(dp) :: ve1(2),ve2(2),ve1p(2),ve2p(2),segdir(2),vl1(2),vl2(2),vint(2),vmerge(2),segnormal(2)
real(dp) :: edge_intersect_norm(cm2dopt%NintEmax)
real(dp) :: edge_intersect(cm2dopt%NintEmax,2)

!Initialise 
sint_idx = 0 
Neintersect = 0 
etype_check = 0 
nselected = 0
zxmin = 0.0d0 
zxmax = 0.0d0 
zymin = 0.0d0 
zymax = 0.0d0 
zzmin = 0.0d0 
zzmax = 0.0d0 
node_select(:) = 0 
edge_check(:) = 0 

!Allocate edge_clip structure 
allocate(edge_clip(volume_mesh_full%nedge))
do ee=1,volume_mesh_full%nedge
    edge_clip(ee)%type = 0 
    edge_clip(ee)%nint = 0
end do 

!Find intersections and classify edges
do ee=1,volume_mesh_full%nedge

    !Edge vertices 
    ev1 = volume_mesh_full%edge(ee,1)
    ev2 = volume_mesh_full%edge(ee,2)

    !If either vertex belongs to a type 1 surface intersecting cell then check this edge 
    echeck = 0 
    if ((vtx_type1(ev1) == 1) .OR. (vtx_type1(ev2) == 1)) then 
        echeck = 1
    end if 

    !Clip edge if required 
    if (echeck == 0) then !Set as fully intetnal or external based on end vertex states 
        if ((vtx_external(ev1) == 1) .AND. (vtx_external(ev2) == 1)) then !set external
            edge_clip(ee)%type = 2
        elseif ((vtx_external(ev1) == 0) .AND. (vtx_external(ev2) == 0)) then !set internal 
            edge_clip(ee)%type = 3
        else !error edge type
            print *, '** unkonwn edge case - no clip check'
            cm2dopt%cm2dfailure = 1
            return 
        end if 
    else !Clip edge

        !Edge end vertices
        ve1(:) = volume_mesh_full%vertices(ev1,:)
        ve2(:) = volume_mesh_full%vertices(ev2,:)

        !Pad segment length 
        segdir(:) = ve2(:) - ve1(:)
        segnorm = norm2(segdir(:))
        ve2p(:) = ve2(:) + segdir(:)*cm2dopt%elenpad
        ve1p(:) = ve1(:) - segdir(:)*cm2dopt%elenpad

        !Padding size 
        cpadSZ = real(segnorm,dp)

        !Intersection bounding box
        zxmin = min(ve1p(1),ve2p(1)) - cpadSZ*cm2dopt%ADTpadding !tgt bounding box -> xmin
        zxmax = max(ve1p(1),ve2p(1)) + cpadSZ*cm2dopt%ADTpadding !tgt bounding box -> xmax
        zymin = min(ve1p(2),ve2p(2)) - cpadSZ*cm2dopt%ADTpadding !tgt bounding box -> ymin
        zymax = max(ve1p(2),ve2p(2)) + cpadSZ*cm2dopt%ADTpadding !tgt bounding box -> ymax

        !Identify any edge bounding boxes that may overlap the cell 
        call search_ADtree(nselected,node_select,surface_adtree,zxmin,zxmax,zymin,zymax,zzmin,zzmax)

        !If any nodes are identified then search these for intersections 
        N_edge_intersect = 0 
        if (nselected .NE. 0) then 
            edge_intersect_seg(:) = 0 
            edge_intersect_type(:) = 0 
            edge_intersect(:,:) = 0.0d0 
            do nn=1,nselected
                do kk=1,surface_adtree%tree(node_select(nn))%nentry
    
                    !Verticies of the ends of the surface segment 
                    vl1(:) = surface_mesh%vertices(surface_mesh%faces(surface_adtree%tree(node_select(nn))%entry(kk),1),:)
                    vl2(:) = surface_mesh%vertices(surface_mesh%faces(surface_adtree%tree(node_select(nn))%entry(kk),2),:)

                    !Intersection type if any 
                    segval = seg_seg_intersect_bool(vl1,vl2,ve1,ve2)
                    if (segval .NE. 0) then 
                    ! if (segval == 1) then 

                        !Find intersection location (is within edge and surface segment here)
                        vint = line_line_intersection_loc_inl1(ve1,ve2,vl1,vl2)

                        !Test against existing intersections to ensure it is new and valid --- 
                        !Check if a new surface segment
                        segnew = 1 
                        do ii=1,N_edge_intersect
                            if (surface_adtree%tree(node_select(nn))%entry(kk) == edge_intersect_seg(N_edge_intersect)) then 
                                segnew = 0 
                            end if
                        end do  

                        !Check if co-incidint with another intersection on this edge 
                        intnew = 1
                        do ii=1,N_edge_intersect
                            if (norm2(vint(:)-edge_intersect(ii,:)) .LE. set_cmp_prc(segnorm*cm2dopt%intcointol,min_precision)) then
                                intnew = 0 
                            end if 
                        end do 

                        !If near an edge end then check if co-incidint with an intersection on another edge 
                        int_unique = 0 
                        efrac = norm2(vint(:) -  ve1(:))/segnorm
                        if ((efrac .GT. 0.99d0) .OR. (efrac .LT. 0.01d0)) then 
                            if (efrac .GT. 0.5d0) then 
                                vref = ev2
                            else
                                vref = ev1
                            end if 
                            do aa=1,4 
                                eadj = volume_mesh_full%V2E(vref,aa)
                                if ((eadj .GT. 0) .AND. (eadj .NE. ee)) then 
                                    do ii=1,edge_clip(eadj)%nint
                                        idist = norm2(vint(:) - edge_clip(eadj)%intloc(ii,:))
                                        if (idist .LE. set_cmp_prc(segnorm*cm2dopt%intcointol,min_precision)) then 
                                            int_unique = edge_clip(eadj)%vtx_idx(ii)
                                            vmerge(:) = edge_clip(eadj)%intloc(ii,:)
                                            exit 
                                        end if
                                    end do 
                                    if (int_unique .GT. 0) then 
                                        exit 
                                    end if 
                                end if 
                            end do 
                        end if 

                        !Add intersection if valid and new 
                        if ((segnew == 1) .AND. (intnew == 1)) then 
                            N_edge_intersect = N_edge_intersect + 1
                            if (N_edge_intersect .GT. cm2dopt%NintEmax) then 
                                cm2dopt%cm2dfailure = 1
                                print *, '** maximum number of edge-geometry intersections exceeded on edge ',ee
                                return 
                            end if 
                            if (int_unique == 0) then !Unique or new so add
                                edge_intersect_type(N_edge_intersect) = -1
                                edge_intersect(N_edge_intersect,:) = vint(:)
                                edge_intersect_seg(N_edge_intersect) = surface_adtree%tree(node_select(nn))%entry(kk)
                            else !Merge with vertex int_unique
                                edge_intersect_type(N_edge_intersect) = int_unique
                                edge_intersect(N_edge_intersect,:) = vmerge(:)
                                edge_intersect_seg(N_edge_intersect) = surface_adtree%tree(node_select(nn))%entry(kk)
                            end if 
                        end if
                    end if 
                end do 
            end do 
        end if 
    
        !If any intersections 
        if (N_edge_intersect .GT. 0) then !Order intersections and identify internal and external sections of the edge 


            !Count intersecting edge 
            Neintersect = Neintersect + 1

            !Set edge type to intersecting 
            edge_clip(ee)%type = 4

            !Store count
            edge_clip(ee)%nint = N_edge_intersect

            !Store reference end 1 of edge 
            edge_clip(ee)%refend1 = ev1

            !Set intersection list to ordered from edge start refend1
            int_selected(:) = 0 
            edge_intersect_norm(:) = 0.0d0  
            do kk=1,N_edge_intersect
                edge_intersect_norm(kk) = norm2(edge_intersect(kk,:) - ve1(:))
            end do 
            allocate(edge_clip(ee)%intloc(N_edge_intersect,2))
            allocate(edge_clip(ee)%vtx_idx(N_edge_intersect))
            allocate(edge_clip(ee)%inttype(N_edge_intersect))
            allocate(edge_clip(ee)%intidx(N_edge_intersect))  
            allocate(edge_clip(ee)%surfseg(N_edge_intersect)) 
            allocate(edge_clip(ee)%intfrac(N_edge_intersect)) 
            edge_clip(ee)%vtx_idx(:) = 0 !construct later when these vertices are added to the mesh
            do kk=1,N_edge_intersect

                !Select intersection at maximum distance
                intselect = maxloc(edge_intersect_norm(1:N_edge_intersect),1,&
                                   mask = int_selected(1:N_edge_intersect) == 0)

                !Store 
                edge_clip(ee)%intloc(N_edge_intersect-kk+1,:) = edge_intersect(intselect,:)
                edge_clip(ee)%inttype(N_edge_intersect-kk+1) = edge_intersect_type(intselect)
                sint_idx = sint_idx + 1
                edge_clip(ee)%intidx(N_edge_intersect-kk+1) = sint_idx
                edge_clip(ee)%surfseg(N_edge_intersect-kk+1) = edge_intersect_seg(intselect)
                edge_clip(ee)%intfrac(N_edge_intersect-kk+1) = real(edge_intersect_norm(intselect)/segnorm,dp) 

                ! !Increment vertex index if not a merge (vertex indecies are not needed here)
                ! if (edge_intersect_type(intselect) .GT. 0) then !merge
                !     edge_clip(ee)%vtx_idx(N_edge_intersect-kk+1) = edge_intersect_type(intselect)
                ! else !normal
                !     ! vtx_idx_smesh = vtx_idx_smesh + 1 
                !     ! edge_clip(ee)%vtx_idx(N_edge_intersect-kk+1) = vtx_idx_smesh
                ! end if 

                !Set distance to zero and set to visited
                int_selected(intselect) = 1
                edge_intersect_norm(intselect) = 0.0d0 
            end do 

            !Classify each intersection as entry or exit 
            allocate(edge_clip(ee)%int_inout(N_edge_intersect))
            edge_clip(ee)%int_inout(:) = 0 
            do kk=1,N_edge_intersect
                vl1(:) = surface_mesh%vertices(surface_mesh%faces(edge_clip(ee)%surfseg(kk),1),:)
                vl2(:) = surface_mesh%vertices(surface_mesh%faces(edge_clip(ee)%surfseg(kk),2),:)
                segnormal(1) = vl1(2) - vl2(2)
                segnormal(2) = vl2(1) - vl1(1)
                if (dot_product(segnormal,segdir) .GT. 0.0d0) then !Exit 
                    edge_clip(ee)%int_inout(kk) = -1
                elseif (dot_product(segnormal,segdir) .LT. 0.0d0) then !Entry 
                    edge_clip(ee)%int_inout(kk) = 1
                else
                    cm2dopt%cm2dfailure = 1
                    print *, '** unable to identify state of edge surface intersection'
                    return 
                end if 
            end do 

            !Classify internal/external status of each segment of this edge ----------
            !Allocate
            allocate(edge_clip(ee)%int_seg_type(N_edge_intersect+1))
            edge_clip(ee)%int_seg_type(:) = 0 

            !Set first segment 
            if ((edge_clip(ee)%int_inout(1) == 1) .AND. (vtx_external(ev1) == 1)) then !entry so segment external 
                edge_clip(ee)%int_seg_type(1) = 2
            elseif ((edge_clip(ee)%int_inout(1) == -1) .AND. (vtx_external(ev1) == 0)) then !exit so segment internal 
                edge_clip(ee)%int_seg_type(1) = 3
            ! elseif (vtx_external(ev1) == 1) then !fallback set as external if first vertex is external
            !     edge_clip(ee)%int_seg_type(1) = 2
            ! elseif (vtx_external(ev1) == 0) then !fallback set as internal if first vertex is internal
            !     edge_clip(ee)%int_seg_type(1) = 3
            else 
                cm2dopt%cm2dfailure = 1
                write(*,'(A,I0)') '    ** edge intersection enter/exit disagreement with vtxinternal at edge: ',ee
                print *, edge_clip(ee)%nint
                print *, edge_clip(ee)%int_inout(1)
                print *, edge_clip(ee)%type
                print *, vtx_external(ev1)
                print *, volume_mesh_full%vertices(ev1,:)
                print *, volume_mesh_full%vertices(ev2,:)
                return
            end if 

            !Set remaining segments 
            do kk=2,N_edge_intersect+1
                if (edge_clip(ee)%int_seg_type(kk-1) == 2) then !currently external
                    if (edge_clip(ee)%int_inout(kk-1) == 1) then !entry -> current segment internal
                        edge_clip(ee)%int_seg_type(kk) = 3
                    elseif (edge_clip(ee)%int_inout(kk-1) == -1) then !exit from external -> warn (keep state)
                        edge_clip(ee)%int_seg_type(kk) = 2
                        print *, '** warning: edge exiting geometry from external state '
                        return 
                        cm2dopt%cm2dfailure = 1
                    elseif (edge_clip(ee)%int_inout(kk-1) == 0) then !no transition so keep state
                        edge_clip(ee)%int_seg_type(kk) = 2
                    end if 
                elseif (edge_clip(ee)%int_seg_type(kk-1) == 3) then !currently internal
                    if (edge_clip(ee)%int_inout(kk-1) == 1) then !entry -> warn (keep state)
                        edge_clip(ee)%int_seg_type(kk) = 3
                        print *, '** warning: edge entering geometry from internal state '
                        cm2dopt%cm2dfailure = 1
                        return 
                    elseif (edge_clip(ee)%int_inout(kk-1) == -1) then !exit from internal ->  current segment external
                        edge_clip(ee)%int_seg_type(kk) = 2
                    elseif (edge_clip(ee)%int_inout(kk-1) == 0) then !no transition so keep state
                        edge_clip(ee)%int_seg_type(kk) = 3
                    end if 
                end if 
            end do 

            !Build edge intersection mesh ----------
            !Allocate
            allocate(edge_clip(ee)%edge_mesh(N_edge_intersect+1,2))
            
            !Set first segment (nagative numbers index the intersections along the edge)
            edge_clip(ee)%edge_mesh(1,1) = ev1
            edge_clip(ee)%edge_mesh(1,2) = -1

            !Set middle segments 
            if (N_edge_intersect .GE. 2) then 
                do kk=2,N_edge_intersect
                    edge_clip(ee)%edge_mesh(kk,1) = -1*(kk - 1)
                    edge_clip(ee)%edge_mesh(kk,2) = -1*kk
                end do 
            end if 

            !Set final segment 
            edge_clip(ee)%edge_mesh(N_edge_intersect+1,1) = -N_edge_intersect
            edge_clip(ee)%edge_mesh(N_edge_intersect+1,2) = ev2
        else !Set as fully intetnal or external based on end vertex states 
            if ((vtx_external(ev1) == 1) .AND. (vtx_external(ev2) == 1)) then !set external
                edge_clip(ee)%type = 2
            elseif ((vtx_external(ev1) == 0) .AND. (vtx_external(ev2) == 0)) then !set internal 
                edge_clip(ee)%type = 3
            else !Set edge type to check subsequently      !error edge type
                edge_clip(ee)%type = 5
                etype_check = etype_check + 1
                edge_check(etype_check) = ee 
            end if 
        end if 
    end if 
end do 

!Check ambigous edges and set vtx_external on these edges if required 
if (etype_check .NE. 0) then 
    call set_amb_edge_state(edge_clip,vtx_external,etype_check,edge_check,volume_mesh_full%edge)
end if 

!Debug ================================
! !Debug write error edges -------
! print *, 'Necheck = ',etype_check
! open(11,file='io/edgeerror.dat')
! do ee=1,etype_check
!     kk = edge_check(ee)
!     ev1 = vmf_edges(kk,1)
!     ev2 = vmf_edges(kk,2)
!     write(11,*) volume_mesh_full%vtx(ev1,:) , volume_mesh_full%vtx(ev2,:)
! end do 
! close(11)

! !Debug write reference external vertices ------- 
! open(11,file='io/vtxexternal.dat')
! do ee=1,volume_mesh_full%nvtx
!     if (vtx_external(ee) == 1) then 
!         write(11,*) volume_mesh_full%vtx(ee,:)
!     end if
! end do 
! close(11)
return 
end subroutine clip_edges2surface




!Check and set state of ambugous edges subroutine ===========================
subroutine set_amb_edge_state(edge_clip,vtx_external,etype_check,edge_check,vmf_edges)
implicit none 

!Variables - Import
integer(in) :: etype_check
integer(in), dimension(:) :: edge_check,vtx_external
integer(in), dimension(:,:) :: vmf_edges
type(edgeint), dimension(:), allocatable :: edge_clip

!Variables - Local 
integer(in) :: ff,ee 
integer(in) :: nupdate,etgt,ev1,ev2

!Flood check ambugous edges 
do ff=1,etype_check !10*etype_check
    nupdate = 0 
    do ee=1,etype_check
        etgt = edge_check(ee)
        if (edge_clip(etgt)%type == 5) then 
            nupdate = nupdate + 1
            ev1 = vmf_edges(etgt,1)
            ev2 = vmf_edges(etgt,2)
            if (vtx_external(ev1) .NE. vtx_external(ev2)) then 
                if ((vtx_external(ev1) == 1) .OR. (vtx_external(ev2) == 1)) then 
                    edge_clip(etgt)%type = 2
                    vtx_external(ev1) = 1
                    vtx_external(ev2) = 1
                end if  
            elseif ((vtx_external(ev1) == 1) .AND. (vtx_external(ev2) == 1)) then 
                edge_clip(etgt)%type = 2
            end if 
        end if 
    end do 

    !Exit if no more edges to update 
    if (nupdate == 0) then 
        exit 
    end if 
end do 
return 
end subroutine set_amb_edge_state




!Build base full mesh subroutine ===========================
subroutine build_full_mesh_EXACT(volume_mesh_full,qt_mesh,cell_keep,cell_2_qtcell,eopposite,celledges)
implicit none 

!Variables - Import
integer(in) :: eopposite(4),celledges(4,2)
integer(in), dimension(:) :: cell_keep,cell_2_qtcell
type(quadtree_data) :: qt_mesh
type(vol_mesh_data) :: volume_mesh_full 

!Variables - Local 
integer(in) :: cc,aa,Nvtx,Nedge,Ncell,c1,c2,eins 
integer(in) :: cadj,cadj_valid,build_edg
integer(in) :: vtx_index(qt_mesh%vins-1),cell_index(qt_mesh%cins-1),cell_edg_idx(qt_mesh%cins-1,4)

!Index mesh edges 
Nvtx = 0 
Nedge = 0 
Ncell = 0 
vtx_index(:) = 0 
cell_index(:) = 0
cell_2_qtcell(:) = 0  
cell_edg_idx(:,:) = 0 
do cc=1,qt_mesh%cins-1
    if (cell_keep(cc) .NE. 0) then 
        Ncell = Ncell + 1
        cell_index(cc) = Ncell
        cell_2_qtcell(Ncell) = cc
        do aa=1,4

            !Adjacent cell 
            cadj = qt_mesh%cell_adjacent(cc,aa)

            !Check valid cell and edge 
            cadj_valid = 0 
            if (cadj .GT. 0) then 
                if (qt_mesh%cell_level(cadj) .LE. qt_mesh%cell_level(cc)) then 
                    if (cell_keep(cadj) .NE. 0) then 
                        cadj_valid = 1
                    end if 
                end if 
            elseif (cadj == -2) then 
                cadj_valid = 1
            end if

            !If valid edge check if new and tag to build if required
            build_edg = 0 
            if (cadj_valid == 1) then 
                if (cell_edg_idx(cc,aa) == 0) then 
                    build_edg = 1
                end if
            end if 
            
            !Build edge if required 
            if (build_edg == 1) then 

                !Check edge end vertices v1/v2 and build if required 
                if (vtx_index(qt_mesh%cell_vcnr(cc,celledges(aa,1))) == 0) then !Build 
                    Nvtx = Nvtx + 1
                    vtx_index(qt_mesh%cell_vcnr(cc,celledges(aa,1))) = Nvtx
                end if
                if (vtx_index(qt_mesh%cell_vcnr(cc,celledges(aa,2))) == 0) then !Build 
                    Nvtx = Nvtx + 1
                    vtx_index(qt_mesh%cell_vcnr(cc,celledges(aa,2))) = Nvtx
                end if

                !Build edge 
                Nedge = Nedge + 1
                cell_edg_idx(cc,aa) = Nedge

                !If across equal refinement level boundary push edge to adjacent cell 
                if (cadj .GT. 0) then 
                    if (qt_mesh%cell_level(cadj) == qt_mesh%cell_level(cc)) then 
                        cell_edg_idx(cadj,eopposite(aa)) = cell_edg_idx(cc,aa)
                    end if
                end if 
            end if 
        end do 
    end if
end do 

!Tag vertices that belong to any geometry external cells 
allocate(volume_mesh_full%vtx_type(Nvtx))
volume_mesh_full%vtx_type(:) = 0 
do cc=1,qt_mesh%cins-1
    if (cell_keep(cc) == 2) then 
        volume_mesh_full%vtx_type(qt_mesh%cell_vcnr(cc,:)) = 1 
    end if  
end do 

!Allocate mesh 
volume_mesh_full%nvtx = Nvtx
allocate(volume_mesh_full%vertices(volume_mesh_full%nvtx,2))
volume_mesh_full%vertices(:,:) = qt_mesh%vtx(1:Nvtx,:) 
allocate(volume_mesh_full%vtx_surfseg(volume_mesh_full%nvtx))
volume_mesh_full%vtx_surfseg(:) = 0
volume_mesh_full%nedge = Nedge
allocate(volume_mesh_full%edge(volume_mesh_full%nedge,4))
volume_mesh_full%edge(:,:) = 0
volume_mesh_full%ncell = Ncell
allocate(volume_mesh_full%cell_level(volume_mesh_full%ncell))
do cc=1,qt_mesh%cins-1
    if (cell_index(cc) .GT. 0) then 
        volume_mesh_full%cell_level(cell_index(cc)) = qt_mesh%cell_level(cc)
    end if
end do 

!Build full mesh 
eins = 0 
do cc=1,qt_mesh%cins-1
    do aa=1,4
        if (cell_edg_idx(cc,aa) .GT. 0) then 

            !Adjacent cell 
            cadj = qt_mesh%cell_adjacent(cc,aa)

            !Cells 
            c1 = cell_index(cc)
            if (cadj .GT. 0) then 
                c2 = cell_index(cadj)
            else
                c2 = cadj
            end if 

            !Build edge 
            eins = eins + 1
            volume_mesh_full%edge(eins,1) = qt_mesh%cell_vcnr(cc,celledges(aa,1))
            volume_mesh_full%edge(eins,2) = qt_mesh%cell_vcnr(cc,celledges(aa,2))
            volume_mesh_full%edge(eins,3) = c2 
            volume_mesh_full%edge(eins,4) = c1 

            !Tag edge as constructed
            cell_edg_idx(cc,aa) = -cell_edg_idx(cc,aa)
            if (cadj .GT. 0) then 
                if (qt_mesh%cell_level(cadj) == qt_mesh%cell_level(cc)) then 
                    cell_edg_idx(cadj,eopposite(aa)) = cell_edg_idx(cc,aa)
                end if
            end if 
        end if 
    end do 
end do 
return 
end subroutine build_full_mesh_EXACT


end module cellmesh2d_mesh_build_exact_mod