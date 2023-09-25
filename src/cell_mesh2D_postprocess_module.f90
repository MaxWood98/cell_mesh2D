!cell_mesh2d mesh postprocessing module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 5.0
!Updated 06-09-2023

!Module
module cellmesh2d_postprocess_mod
use cellmesh2d_geometry_mod
contains  


!Laplacian near surface smoothing subroutine ===========================
subroutine nearsurf_lap_smooth(volume_mesh,cm2dopt,MaxValence)
implicit none 

!Variables - Import
integer(in) :: MaxValence
type(cm2d_options) :: cm2dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: aa,ee,ff,vv,ss,flevel,Ntag,v1,v2  
integer(in) :: vtx_type(volume_mesh%nvtx),vtxidxF(volume_mesh%nvtx)
integer(in), dimension(:), allocatable :: vtxidxR,nvatt
integer(in), dimension(:,:), allocatable :: v2v
real(dp), dimension(:,:), allocatable :: vposN

!Display 
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '--> laplacian smothing near surfaces'
end if

!Tag vertices on object surfaces
vtx_type(:) = 0 
do ee=1,volume_mesh%nedge
    if (volume_mesh%edge(ee,3) == -1) then 
        vtx_type(volume_mesh%edge(ee,1)) = 1
        vtx_type(volume_mesh%edge(ee,2)) = 1
    end if 
end do 

!Flood vertex types away from surface
flevel = 2
do ff=1,cm2dopt%nlpflood
    do ee=1,volume_mesh%nedge
        if (vtx_type(volume_mesh%edge(ee,1)) == flevel-1) then 
            if (vtx_type(volume_mesh%edge(ee,2)) == 0) then 
                vtx_type(volume_mesh%edge(ee,2)) = flevel
            end if
        end if
        if (vtx_type(volume_mesh%edge(ee,2)) == flevel-1) then 
            if (vtx_type(volume_mesh%edge(ee,1)) == 0) then 
                vtx_type(volume_mesh%edge(ee,1)) = flevel
            end if
        end if
    end do
    flevel = flevel + 1 
end do 

!Untag any vertices on the boundaries to set as fixed 
do ee=1,volume_mesh%nedge
    if ((volume_mesh%edge(ee,3) .LT. 0) .OR. (volume_mesh%edge(ee,4) .LT. 0)) then 
        vtx_type(volume_mesh%edge(ee,1:2)) = 0 
    end if 
end do 

!Index flooded vertices of flevel >= 2
Ntag = 0 
do vv=1,volume_mesh%nvtx
    if (vtx_type(vv) .GE. 2) then 
        Ntag = Ntag + 1
    end if
end do 
if (cm2dopt%dispt == 1) then
    write(*,'(A,I0,A)') '    {identified ',Ntag,' near surface vertices for smoothing}'
end if 
allocate(vtxidxR(Ntag))
vtxidxR(:) = 0 
vtxidxF(:) = 0 
Ntag = 0 
do vv=1,volume_mesh%nvtx
    if (vtx_type(vv) .GE. 2) then 
        Ntag = Ntag + 1
        vtxidxR(Ntag) = vv 
        vtxidxF(vv) = Ntag 
    end if
end do

!Build v2v for tagged vertices 
allocate(v2v(Ntag,MaxValence))
allocate(nvatt(Ntag))
v2v(:,:) = 0 
nvatt(:) = 0 
do ee=1,volume_mesh%nedge

    !Edge vertices 
    v1 = volume_mesh%edge(ee,1)
    v2 = volume_mesh%edge(ee,2)

    !Add opposite vertex if one is tagged 
    if (vtxidxF(v1) .NE. 0) then !add v2 to v2v for v1
        do aa=1,MaxValence
            if (v2v(vtxidxF(v1),aa) == 0) then 
                v2v(vtxidxF(v1),aa) = v2
                nvatt(vtxidxF(v1)) = nvatt(vtxidxF(v1)) + 1
                exit 
            end if
        end do 
    end if
    if (vtxidxF(v2) .NE. 0) then !add v1 to v2v for v2
        do aa=1,MaxValence
            if (v2v(vtxidxF(v2),aa) == 0) then 
                v2v(vtxidxF(v2),aa) = v1
                nvatt(vtxidxF(v2)) = nvatt(vtxidxF(v2)) + 1
                exit 
            end if
        end do 
    end if
end do 

!Smooth tagged vertices 
allocate(vposN(Ntag,2))
do ss=1,cm2dopt%nlpsmooth
    do vv=1,Ntag
        vposN(vv,:) = 0.0d0 
        do aa=1,nvatt(vv)
            vposN(vv,1) = vposN(vv,1) + volume_mesh%vertices(v2v(vv,aa),1)
            vposN(vv,2) = vposN(vv,2) + volume_mesh%vertices(v2v(vv,aa),2)
        end do 
        vposN(vv,:) = vposN(vv,:)/real(nvatt(vv),dp)
    end do 
    do vv=1,Ntag
        volume_mesh%vertices(vtxidxR(vv),:) = vposN(vv,:)
    end do 
end do 
return 
end subroutine nearsurf_lap_smooth




!Mesh short edge cleaning subroutine ===========================
subroutine clean_mesh_shortE(volume_mesh,cm2dopt)
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ii,vv
integer(in) :: v1,v2,cl,cr,lr,ll,Neshort,Lref,NvmergeN,Nvkeep,Nvnew,Nenew,Nvtx0,Nset_merge,Nmset_final,stgt
integer(in) :: edge_status(volume_mesh%nedge),vtx_status(volume_mesh%nvtx),set_merge(volume_mesh%nvtx)
integer(in) :: vtx_map(volume_mesh%nvtx),merge_set_map(volume_mesh%nvtx)
integer(in) :: edge_new(volume_mesh%nedge,4)
integer(in) :: vtx_surfsegN(volume_mesh%nvtx)
real(dp) :: edgelen,edgelen_ref
real(dp) :: vtx_new(volume_mesh%nvtx,2)

!Identify any short edges 
Neshort = 0 
edge_status(:) = 0 !0 = normal | 1 = short
do ii=1,volume_mesh%nedge
    
    !Vertices on edge
    v1 = volume_mesh%edge(ii,1)
    v2 = volume_mesh%edge(ii,2)

    !Cells on edge 
    cl = volume_mesh%edge(ii,3)
    cr = volume_mesh%edge(ii,4)

    !Find maximum adjacent level 
    lr = volume_mesh%cell_level(cr)
    if (cl .GT. 0) then 
        ll = volume_mesh%cell_level(cl)
    else
        ll = 0 
    end if
    Lref = max(lr,ll)

    !Find reference edge length for this cell 
    edgelen_ref = cm2dopt%EminLength*(2.0d0*cm2dopt%far_field_bound/(2.0d0**(Lref - 1)))

    !Test if edge is too short
    edgelen = norm2(volume_mesh%vertices(v2,:) - volume_mesh%vertices(v1,:))
    if (edgelen .LE. edgelen_ref) then 
        Neshort = Neshort + 1
        edge_status(ii) = 1
    end if
end do 

!Clean if any short edges are identified 
if (Neshort .NE. 0) then 

    !Properties to store
    Nvtx0 = volume_mesh%nvtx

    !Find vertices to merge 
    NvmergeN = 0 
    Nset_merge = 0 
    set_merge(:) = 0 
    vtx_status(:) = 0 
    do ii=1,volume_mesh%nedge
        if (edge_status(ii) == 1) then 

            !Vertices on edge
            v1 = volume_mesh%edge(ii,1)
            v2 = volume_mesh%edge(ii,2)

            !Set merge status 
            if ((vtx_status(v1) == 0) .AND. (vtx_status(v2) == 0)) then 
                NvmergeN = NvmergeN + 1
                vtx_status(v1) = NvmergeN
                vtx_status(v2) = NvmergeN
            elseif ((vtx_status(v1) == 0) .AND. (vtx_status(v2) .NE. 0)) then 
                vtx_status(v1) = vtx_status(v2)
            elseif ((vtx_status(v1) .NE. 0) .AND. (vtx_status(v2) == 0)) then 
                vtx_status(v2) = vtx_status(v1)
            else !Both are non zero merge sets
                if (vtx_status(v1) == vtx_status(v2)) then 
                    !Do nothing as the two vertices are already in the same merge set
                else !Must merge two full merge sets 
                    if ((set_merge(vtx_status(v1)) .NE. 0) .AND. &
                        (set_merge(vtx_status(v2)) == 0)) then !merge into existing set
                        set_merge(vtx_status(v2)) = set_merge(vtx_status(v1))
                    elseif ((set_merge(vtx_status(v1)) == 0) .AND. &
                            (set_merge(vtx_status(v2)) .NE. 0)) then !merge into existing set
                        set_merge(vtx_status(v1)) = set_merge(vtx_status(v2))
                    elseif ((set_merge(vtx_status(v1)) == 0) .AND. &
                            (set_merge(vtx_status(v2)) == 0)) then !neither set has been merge tagged so merge both into set 1 
                        Nset_merge = Nset_merge + 1
                        ! set_merge(vmerge(v1)) = Nset_merge
                        ! set_merge(vmerge(v2)) = Nset_merge
                        set_merge(vtx_status(v1)) = vtx_status(v1)
                        set_merge(vtx_status(v2)) = vtx_status(v1)
                    else 
                        if (set_merge(vtx_status(v1)) == set_merge(vtx_status(v2))) then 
                            !Do nothing as both sets are already tagged to merge 
                        else
                            print *, '** vertex merging error on short edge merging : ',&
                            vtx_status(v1),' // ',vtx_status(v2)
                            print *,set_merge(vtx_status(v1)),' // ',set_merge(vtx_status(v2))
                        end if 
                    end if 
                end if
            end if
        end if
    end do 

    !Merge all merge sets to build final set of merge sets 
    Nmset_final = 0 
    merge_set_map(:) = 0 
    do vv=1,volume_mesh%nvtx
        if (vtx_status(vv) .NE. 0) then 
            stgt = vtx_status(vv)
            if (set_merge(stgt) .NE. 0) then 
                if ((merge_set_map(set_merge(stgt)) == 0) .AND. (merge_set_map(stgt) == 0)) then 
                    Nmset_final = Nmset_final + 1
                    merge_set_map(stgt) = Nmset_final
                    merge_set_map(set_merge(stgt)) = Nmset_final
                elseif ((merge_set_map(set_merge(stgt)) .NE. 0) .AND. (merge_set_map(stgt) == 0)) then 
                    merge_set_map(stgt) = merge_set_map(set_merge(stgt))
                elseif ((merge_set_map(set_merge(stgt)) == 0) .AND. (merge_set_map(stgt) .NE. 0)) then  
                    merge_set_map(set_merge(stgt)) = merge_set_map(stgt)
                end if 
            elseif (merge_set_map(stgt) == 0) then 
                Nmset_final = Nmset_final + 1
                merge_set_map(stgt) = Nmset_final
            end if 
        end if 
    end do

    !Update list of vertex merge sets 
    do vv=1,volume_mesh%nvtx
        if (vtx_status(vv) .NE. 0) then 
            vtx_status(vv) = merge_set_map(vtx_status(vv))
        end if 
    end do 

    !Build new vertex set  
    Nvkeep = 0  
    vtx_map(:) = 0 
    do ii=1,volume_mesh%nvtx 
        if (vtx_status(ii) == 0) then 
            Nvkeep = Nvkeep + 1
            vtx_map(ii) = Nvkeep
        end if
    end do 
    do ii=1,volume_mesh%nvtx 
        if (vtx_status(ii) .NE. 0) then 
            vtx_map(ii) = Nvkeep + vtx_status(ii)
        end if
    end do 
    do ii=1,volume_mesh%nvtx 
        vtx_new(vtx_map(ii),:) = volume_mesh%vertices(ii,:)
        vtx_surfsegN(vtx_map(ii)) = volume_mesh%vtx_surfseg(ii)
    end do 
    Nvnew = Nvkeep + NvmergeN
    volume_mesh%nvtx = Nvnew
    deallocate(volume_mesh%vertices)
    deallocate(volume_mesh%vtx_surfseg)
    allocate(volume_mesh%vertices(volume_mesh%nvtx,2))
    allocate(volume_mesh%vtx_surfseg(volume_mesh%nvtx))
    volume_mesh%vertices(:,:) = vtx_new(1:Nvnew,:)
    volume_mesh%vtx_surfseg(:) = vtx_surfsegN(1:Nvnew)

    !Build new edge set
    Nenew = 0 
    do ii=1,volume_mesh%nedge
        if (edge_status(ii) == 0) then 
            Nenew = Nenew + 1 
            edge_new(Nenew,:) = volume_mesh%edge(ii,:)
            edge_new(Nenew,1) = vtx_map(edge_new(Nenew,1))
            edge_new(Nenew,2) = vtx_map(edge_new(Nenew,2))
        end if 
    end do 
    deallocate(volume_mesh%edge)
    volume_mesh%nedge = Nenew
    allocate(volume_mesh%edge(volume_mesh%nedge,4)) 
    volume_mesh%edge(:,:) = edge_new(1:Nenew,:) 

    !Display
    if (cm2dopt%dispt == 1) then
        write(*,'(A,I0,A)') '    {identified and collapsed ',Neshort,' short edges}'
    end if 
end if
return 
end subroutine clean_mesh_shortE




!Mesh sliver cell cleaning subroutine ===========================
subroutine clean_mesh_sliverC(volume_mesh,cm2dopt,Nmerge,Nmerge_fail)
implicit none 

!Variables - Import
integer(in) :: Nmerge,Nmerge_fail
type(cm2d_options) :: cm2dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ii,aa,aa2,aa3,Ncremove,maxcedge,etgt,ctgt,ctgtp,cadj,cadj2,NcellN,NedgeN,merge_invalid
integer(in) :: Ninvalid,invalid_trigger
integer(in) :: cell_remove(volume_mesh%ncell),cell_surfadj(volume_mesh%ncell),cell_map(volume_mesh%ncell)
integer(in) :: cell_nedge(volume_mesh%ncell),edge_remove(volume_mesh%nedge),edge_map(volume_mesh%nedge)
integer(in), dimension(:), allocatable :: cell_level_new
integer(in), dimension(:,:), allocatable :: cell2edge,cell2cell,edges_new
real(dp) :: AEdge,vref,vmax
real(dp) :: Cvol(volume_mesh%ncell)

!Initialise
Nmerge = 0 
Ninvalid = 0 
Nmerge_fail = 0 

!Evaluate cell volumes and tag surface adjacent cells 
Cvol(:) = 0.0d0
cell_surfadj(:) = 0 
do ii=1,volume_mesh%nedge

    !Accumulate volume 
    AEdge = Asegment(volume_mesh%vertices(volume_mesh%edge(ii,1),:),volume_mesh%vertices(volume_mesh%edge(ii,2),:))
    Cvol(volume_mesh%edge(ii,4)) = Cvol(volume_mesh%edge(ii,4)) - AEdge
    if (volume_mesh%edge(ii,3) .GT. 0) then
        Cvol(volume_mesh%edge(ii,3)) = Cvol(volume_mesh%edge(ii,3)) + AEdge
    end if

    !Tag if surface adjacent 
    if (volume_mesh%edge(ii,3) == -1) then
        cell_surfadj(volume_mesh%edge(ii,4)) = 1
    end if
end do

!Check for sliver cells adjacent to the object surface
Ncremove = 0 
cell_remove(:) = 0 
do ii=1,volume_mesh%ncell 
    if (cell_surfadj(ii) == 1) then
        vref = (2.0d0*cm2dopt%far_field_bound/(2.0d0**(volume_mesh%cell_level(ii) - 1)))**2
        if (Cvol(ii) .LE. cm2dopt%CminVol*vref) then 
            Ncremove = Ncremove + 1
            cell_remove(ii) = 1
        end if
    end if 
end do 

!Remove cells if required
if (Ncremove .NE. 0) then 

    !Build cell2edge and cell2cell (of internal edges)
    cell_nedge(:) = 0 
    do ii=1,volume_mesh%nedge
        if ((volume_mesh%edge(ii,3) .GT. 0) .AND. (volume_mesh%edge(ii,4) .GT. 0)) then
            cell_nedge(volume_mesh%edge(ii,4)) = cell_nedge(volume_mesh%edge(ii,4)) + 1
            cell_nedge(volume_mesh%edge(ii,3)) = cell_nedge(volume_mesh%edge(ii,3)) + 1
        end if
    end do 
    maxcedge = maxval(cell_nedge)
    allocate(cell2edge(volume_mesh%ncell,maxcedge))
    allocate(cell2cell(volume_mesh%ncell,maxcedge))
    cell2edge(:,:) = 0 
    cell2cell(:,:) = 0
    cell_nedge(:) = 0 
    do ii=1,volume_mesh%nedge
        if ((volume_mesh%edge(ii,3) .GT. 0) .AND. (volume_mesh%edge(ii,4) .GT. 0)) then
            cell_nedge(volume_mesh%edge(ii,3)) = cell_nedge(volume_mesh%edge(ii,3)) + 1
            cell_nedge(volume_mesh%edge(ii,4)) = cell_nedge(volume_mesh%edge(ii,4)) + 1
            cell2edge(volume_mesh%edge(ii,3),cell_nedge(volume_mesh%edge(ii,3))) = ii 
            cell2edge(volume_mesh%edge(ii,4),cell_nedge(volume_mesh%edge(ii,4))) = ii 
            cell2cell(volume_mesh%edge(ii,3),cell_nedge(volume_mesh%edge(ii,3))) = volume_mesh%edge(ii,4)
            cell2cell(volume_mesh%edge(ii,4),cell_nedge(volume_mesh%edge(ii,4))) = volume_mesh%edge(ii,3)
        end if 
    end do 

    !Build cell mapping for retained cells
    NcellN = 0 
    cell_map(:) = 0 
    do ii=1,volume_mesh%ncell 
        if (cell_remove(ii) .NE. 1) then 
            NcellN = NcellN + 1
            cell_map(ii) = NcellN
        end if
    end do 

    !Tag edges to remove to eliminate sliver cells and map removed cells to their replacements
    edge_remove(:) = 0 
    do ii=1,volume_mesh%ncell 
        if (cell_remove(ii) == 1) then 

            !Find adjacent cell with largest volume to merge with 
            etgt = 0
            ctgt = 0 
            vmax = 0.0d0 
            invalid_trigger = 0 
            do aa=1,cell_nedge(ii)
                if (Cvol(cell2cell(ii,aa)) .GT. vmax) then 
                    if (cell_remove(cell2cell(ii,aa)) == 0) then 

                        !Targeted cell 
                        ctgtp = cell2cell(ii,aa)

                        !Check if there is a cell adjacent to the target merge cell that also borders the original cell
                        merge_invalid = 0 
                        do aa2=1,cell_nedge(ctgtp) !cells adjacent to the target merge cell 
                            cadj = cell2cell(ctgtp,aa2)
                            if (cadj .NE. ii) then 
                                do aa3=1,cell_nedge(ii) !cells adjacent to the base cell 
                                    cadj2 = cell2cell(ii,aa3)
                                    if (cadj2 .NE. ii) then 
                                        if (cadj2 == cadj) then 
                                            merge_invalid = 1
                                            invalid_trigger = 1
                                            exit 
                                        end if 
                                    end if 
                                end do 
                                if (merge_invalid == 1) then 
                                    exit 
                                end if 
                            end if 
                        end do 

                        !If valid merge then select cell ctgtp
                        if (merge_invalid == 0) then 
                            vmax = Cvol(cell2cell(ii,aa))
                            etgt = cell2edge(ii,aa)    
                            ctgt = cell2cell(ii,aa)
                        end if 
                    end if
                end if
            end do 

            !Count invalid cells 
            if ((ctgt == 0) .AND. (invalid_trigger == 1)) then 
                Ninvalid = Ninvalid + 1
            end if 

            !Tag if valid merge identified
            if (ctgt .GT. 0) then !Valid merge found -> Map cell for replacement and tag edges joining ii and ctgt for removal
                cell_map(ii) = cell_map(ctgt) 
                Nmerge = Nmerge + 1
                edge_remove(etgt) = 1
                do aa=1,cell_nedge(ii)
                    if (cell2cell(ii,aa) == ctgt) then 
                        edge_remove(cell2edge(ii,aa)) = 1
                    end if 
                end do
            else !No valid merge case identified
                cell_remove(ii) = 2
                Nmerge_fail = Nmerge_fail + 1
            end if 
        end if 
    end do 

    !Update mapping for retained cells 
    do ii=1,volume_mesh%ncell 
        if (cell_remove(ii) == 2) then 
            NcellN = NcellN + 1
            cell_map(ii) = NcellN 
        end if
    end do

    !Tag edges that map to the same cell on both sides to be removed
    do ii=1,volume_mesh%nedge
        if ((volume_mesh%edge(ii,3) .GT. 0) .AND. (volume_mesh%edge(ii,4) .GT. 0)) then 
            if (cell_map(volume_mesh%edge(ii,3)) == cell_map(volume_mesh%edge(ii,4))) then 
                edge_remove(ii) = 1
            end if
        end if
    end do 

    !Build edge mapping for retained edges
    NedgeN = 0 
    edge_map(:) = 0 
    do ii=1,volume_mesh%nedge
        if (edge_remove(ii) .NE. 1) then 
            NedgeN = NedgeN + 1
            edge_map(ii) = NedgeN
        end if
    end do 

    !Rebuild mesh 
    allocate(edges_new(NedgeN,4))
    allocate(cell_level_new(NcellN))
    do ii=1,volume_mesh%nedge
        if (edge_map(ii) .NE. 0) then 
            edges_new(edge_map(ii),:) = volume_mesh%edge(ii,:)
            edges_new(edge_map(ii),4) = cell_map(edges_new(edge_map(ii),4))
            if (edges_new(edge_map(ii),3) .GT. 0) then 
                edges_new(edge_map(ii),3) = cell_map(edges_new(edge_map(ii),3))
            end if 
        end if
    end do 
    do ii=1,volume_mesh%ncell 
        if (cell_map(ii) .NE. 0) then 
            cell_level_new(cell_map(ii)) = volume_mesh%cell_level(ii)
        end if 
    end do 
    deallocate(volume_mesh%edge)
    allocate(volume_mesh%edge(NedgeN,4))
    deallocate(volume_mesh%cell_level)
    allocate(volume_mesh%cell_level(NcellN))
    volume_mesh%cell_level(:) = cell_level_new(:)
    volume_mesh%edge(:,:) = edges_new(:,:)
    volume_mesh%nedge = NedgeN
    volume_mesh%ncell = NcellN

    !Display
    if (cm2dopt%dispt == 1) then
        write(*,'(A,I0,A,I0,A)') '    {identified ',Ncremove,' sliver cells and merged ',Nmerge,' cells}'
        if (Ninvalid .NE. 0) then 
            write(*,'(A,I0,A)') '    {',Ninvalid,' merges prevented due to double adjacency}'
        end if 
    end if 
end if
return
end subroutine clean_mesh_sliverC




!Mesh bisected cell correction subroutine ===========================
subroutine correct_bisected_cells(volume_mesh,cm2dopt)
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ii,jj,cc,ee,aa,vv,v1,v2,c3,c4,vtgt,vbase,etgt,evalid,exist,maxcvtx,Nbisect,Ncell
integer(in) :: cell_nvtx(volume_mesh%ncell)
integer(in) :: v2v(volume_mesh%nvtx,4),v2e(volume_mesh%nvtx,4)
integer(in) :: cell_tag(volume_mesh%nvtx)
integer(in) :: cell_levelN(2*volume_mesh%ncell)
integer(in), dimension(:,:), allocatable :: cell2vtx 

!Initialise
Nbisect = 0 

!Assign vertices to cells 
cell_nvtx(:) = 0 
do ee=1,volume_mesh%nedge
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    if (c3 .GT. 0) then 
        cell_nvtx(c3) = cell_nvtx(c3) + 1 !Add one as each vertex is seen twice per cell when looking at edges
    end if
    if (c4 .GT. 0) then 
        cell_nvtx(c4) = cell_nvtx(c4) + 1 !Add one as each vertex is seen twice per cell when looking at edges
    end if
end do 
maxcvtx = maxval(cell_nvtx)
cell_nvtx(:) = 0 
allocate(cell2vtx(volume_mesh%ncell,maxcvtx))
cell2vtx(:,:) = 0 
do ee=1,volume_mesh%nedge
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    do aa=1,2
        vtgt = volume_mesh%edge(ee,aa)
        if (c3 .GT. 0) then
            exist = 0 
            do vv=1,maxcvtx
                if (cell2vtx(c3,vv) == vtgt) then 
                    exist = 1
                    exit 
                end if
            end do 
            if (exist == 0) then 
                do vv=1,maxcvtx
                    if (cell2vtx(c3,vv) == 0) then 
                        cell2vtx(c3,vv) = vtgt
                        cell_nvtx(c3) = cell_nvtx(c3) + 1
                        exit 
                    end if
                end do 
            end if 
        end if
        if (c4 .GT. 0) then
            exist = 0 
            do vv=1,maxcvtx
                if (cell2vtx(c4,vv) == vtgt) then 
                    exist = 1
                    exit 
                end if
            end do 
            if (exist == 0) then 
                do vv=1,maxcvtx
                    if (cell2vtx(c4,vv) == 0) then 
                        cell2vtx(c4,vv) = vtgt
                        cell_nvtx(c4) = cell_nvtx(c4) + 1
                        exit 
                    end if
                end do 
            end if 
        end if
    end do 
end do 

!Build v2v and v2e 
v2v(:,:) = 0
v2e(:,:) = 0
do ee=1,volume_mesh%nedge

    !Edge ends
    v1 = volume_mesh%edge(ee,1)
    v2 = volume_mesh%edge(ee,2)

    !Add to v2e 
    !v1
    exist = 0 
    do vv=1,4
        if (v2e(v1,vv) == ee) then 
            exist = 1
            exit 
        end if
    end do 
    if (exist == 0) then 
        do vv=1,4
            if (v2e(v1,vv) == 0) then 
                v2e(v1,vv) = ee
                exit 
            end if
        end do 
    end if 

    !v2
    exist = 0 
    do vv=1,4
        if (v2e(v2,vv) == ee) then 
            exist = 1
            exit 
        end if
    end do 
    if (exist == 0) then 
        do vv=1,4
            if (v2e(v2,vv) == 0) then 
                v2e(v2,vv) = ee
                exit 
            end if
        end do 
    end if 

    !Add to v2v
    !v1
    exist = 0 
    do vv=1,4
        if (v2v(v1,vv) == v2) then 
            exist = 1
            exit 
        end if
    end do 
    if (exist == 0) then 
        do vv=1,4
            if (v2v(v1,vv) == 0) then 
                v2v(v1,vv) = v2
                exit 
            end if
        end do 
    end if 

    !v2
    exist = 0 
    do vv=1,4
        if (v2v(v2,vv) == v1) then 
            exist = 1
            exit 
        end if
    end do 
    if (exist == 0) then 
        do vv=1,4
            if (v2v(v2,vv) == 0) then 
                v2v(v2,vv) = v1
                exit 
            end if
        end do 
    end if 
end do 

!Check for and correct bisected cells 
cell_tag(:) = 0 
cell_levelN(:) = 0 
cell_levelN(1:volume_mesh%ncell) = volume_mesh%cell_level(:)
deallocate(volume_mesh%cell_level)
Ncell = volume_mesh%ncell
do cc=1,volume_mesh%ncell 
    if (cell2vtx(cc,1) .GT. 0) then 

        !Flood all vertices in this cell from the first 
        cell_tag(cell2vtx(cc,1)) = 1
        do ii=1,cell_nvtx(cc) !Flood iteration
            do jj=1,cell_nvtx(cc) !Loop all vertices
                if (cell_tag(cell2vtx(cc,jj)) == 1) then !Flood tag to adjacent vertices in this cell along edges that border this cell 
                    vbase = cell2vtx(cc,jj)
                    do aa=1,4
                        vtgt = v2v(vbase,aa) 
                        etgt = v2e(vbase,aa) 
                        if (vtgt .GT. 0) then 

                            !Check if attached vertex is within cell
                            exist = 0 
                            do ee=1,cell_nvtx(cc)
                                if (vtgt == cell2vtx(cc,ee)) then 
                                    exist = 1
                                end if
                            end do 

                            !Check if edge borders this cell 
                            evalid = 0 
                            if ((volume_mesh%edge(etgt,3) == cc) .OR. (volume_mesh%edge(etgt,4) == cc)) then 
                                evalid = 1
                            end if

                            !Flood tag if within cell
                            if ((exist == 1) .AND. (evalid == 1)) then 
                                cell_tag(vtgt) = 1
                            end if 
                        end if 
                    end do 
                end if 
            end do 
        end do 

        !Check if all vertices in the cell have been tagged
        if (sum(cell_tag(cell2vtx(cc,1:cell_nvtx(cc)))) .NE. cell_nvtx(cc)) then !Cell is bisected

            !Incrment count of bisceted cells 
            Nbisect = Nbisect + 1

            !Increment count of mesh cells 
            Ncell = Ncell + 1

            !Assign level of this new cell 
            cell_levelN(Ncell) = cell_levelN(cc)
        
            !Tag adjacency for this cell of all edges on untagged cell vertices in this cell with the new cell index 
            do ii=1,cell_nvtx(cc) !Flood iteration
                if (cell_tag(cell2vtx(cc,ii)) == 0) then 
                    vbase = cell2vtx(cc,ii)
                    do aa=1,4
                        vtgt = v2v(vbase,aa) 
                        etgt = v2e(vbase,aa) 
                        if (vtgt .GT. 0) then 

                            !Check if attached vertex is within cell
                            exist = 0 
                            do ee=1,cell_nvtx(cc)
                                if (vtgt == cell2vtx(cc,ee)) then 
                                    exist = 1
                                end if
                            end do 

                            !Check if edge borders this cell 
                            evalid = 0 
                            if ((volume_mesh%edge(etgt,3) == cc) .OR. (volume_mesh%edge(etgt,4) == cc)) then 
                                evalid = 1
                            end if

                            !Process if valid
                            if ((exist == 1) .AND. (evalid == 1) .AND. (cell_tag(vtgt) == 0)) then 
                                if (volume_mesh%edge(etgt,3) == cc) then 
                                    volume_mesh%edge(etgt,3) = Ncell
                                end if  
                                if (volume_mesh%edge(etgt,4) == cc) then 
                                    volume_mesh%edge(etgt,4) = Ncell
                                end if  
                            end if 
                        end if
                    end do 
                end if 
            end do 
        end if

        !Reset cell tags
        cell_tag(cell2vtx(cc,1:cell_nvtx(cc))) = 0
    end if 
end do 

!Assign new cell count and cell levels 
volume_mesh%ncell = Ncell
allocate(volume_mesh%cell_level(Ncell))
volume_mesh%cell_level(:) = cell_levelN(1:Ncell)

!Display
if (Nbisect .GT. 0) then 
    if (cm2dopt%dispt == 1) then
        write(*,'(A,I0,A)') '    {identified and split ',Nbisect,' bisected cells}'
    end if 
end if 
return 
end subroutine correct_bisected_cells




!Cell index remapping subroutine ===========================
subroutine remap_cell_indecies(volume_mesh)
implicit none 

!Variables - Import
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: cc,ee,c3,c4,NcellN
integer(in) :: cell_indexN(volume_mesh%ncell)
integer(in), dimension(:), allocatable :: cell_levelN

!Remap cells
NcellN = 0 
cell_indexN(:) = 0 
do ee=1,volume_mesh%nedge
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    if (c3 .GT. 0) then 
        if (cell_indexN(c3) == 0) then
            NcellN = NcellN + 1
            cell_indexN(c3) = NcellN
        end if 
    end if
    if (c4 .GT. 0) then 
        if (cell_indexN(c4) == 0) then
            NcellN = NcellN + 1
            cell_indexN(c4) = NcellN
        end if 
    end if
end do 

!Reasign indecies on edges
do ee=1,volume_mesh%nedge
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    if (c3 .GT. 0) then 
        if (cell_indexN(c3) .GT. 0) then
            volume_mesh%edge(ee,3) = cell_indexN(c3)
        end if 
    end if
    if (c4 .GT. 0) then 
        if (cell_indexN(c4) .GT. 0) then
            volume_mesh%edge(ee,4) = cell_indexN(c4)
        end if 
    end if
end do 

!Reassign cell levels 
allocate(cell_levelN(NcellN))
do cc=1,volume_mesh%ncell
    if (cell_indexN(cc) .GT. 0) then 
        cell_levelN(cell_indexN(cc)) = volume_mesh%cell_level(cc)
    end if 
end do 
deallocate(volume_mesh%cell_level)
allocate(volume_mesh%cell_level(NcellN))
volume_mesh%cell_level(:) = cell_levelN(:)

!Assign new cell count 
volume_mesh%ncell = NcellN
return 
end subroutine remap_cell_indecies




!Subroutine to build vertex to geometry surface links ===========================
subroutine build_surface_links(volume_mesh,surface_mesh)
implicit none 

!Variables - Import
type(vol_mesh_data) :: volume_mesh
type(surface_data) :: surface_mesh

!Variables - Local 
integer(in) :: vv,Nvsurf,segtgt
real(dp) :: segfrac
real(dp) :: vs1(2),vs2(2)

!Count surface linked vertices
Nvsurf = 0 
do vv=1,volume_mesh%nvtx
    if (volume_mesh%vtx_surfseg(vv) .NE. 0) then 
        Nvsurf = Nvsurf + 1
    end if 
end do 

!Allocate arrays
volume_mesh%nvtx_surf = Nvsurf
allocate(volume_mesh%surf_vtx(Nvsurf))
allocate(volume_mesh%surf_vtx_seg(Nvsurf))
allocate(volume_mesh%surf_vtx_segfrac(Nvsurf))
Nvsurf = 0 

!Build arrays
do vv=1,volume_mesh%nvtx
    if (volume_mesh%vtx_surfseg(vv) .NE. 0) then 

        !Increment 
        Nvsurf = Nvsurf + 1

        !Store vertex and segment
        volume_mesh%surf_vtx(Nvsurf) = vv 
        volume_mesh%surf_vtx_seg(Nvsurf) = volume_mesh%vtx_surfseg(vv) 

        !Find segment fraction  
        segtgt = volume_mesh%vtx_surfseg(vv)
        vs1(:) = surface_mesh%vertices(surface_mesh%faces(segtgt,1),:)
        vs2(:) = surface_mesh%vertices(surface_mesh%faces(segtgt,2),:)
        segfrac = norm2(volume_mesh%vertices(vv,:) - vs1(:))/norm2(vs2(:) - vs1(:))
        volume_mesh%surf_vtx_segfrac(Nvsurf) = segfrac
    end if
end do 
return 
end subroutine build_surface_links




!Subroutine to set custom boundary conditions where requested ===========================
subroutine set_custom_bcs(volume_mesh,cm2dopt)
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ii,ee
real(dp) :: xmin,xmax,ymin,ymax

!Each zone on each edge 
do ee=1,volume_mesh%nedge
    if ((volume_mesh%edge(ee,3) == -2) .OR. (volume_mesh%edge(ee,4) == -2)) then 
        xmin = minval(volume_mesh%vertices(volume_mesh%edge(ee,1:2),1))
        xmax = maxval(volume_mesh%vertices(volume_mesh%edge(ee,1:2),1))
        ymin = minval(volume_mesh%vertices(volume_mesh%edge(ee,1:2),2))
        ymax = maxval(volume_mesh%vertices(volume_mesh%edge(ee,1:2),2))
        do ii=1,cm2dopt%Nzone_cBC
            if ((xmin .LE. cm2dopt%BC_zone_coords(ii,2)) .AND. (xmax .GE. cm2dopt%BC_zone_coords(ii,1))) then 
                if ((ymin .LE. cm2dopt%BC_zone_coords(ii,4)) .AND. (ymax .GE. cm2dopt%BC_zone_coords(ii,3))) then 
                    if (volume_mesh%edge(ee,3) == -2) then 
                        volume_mesh%edge(ee,3) = cm2dopt%BC_zone_bc(ii)
                    elseif (volume_mesh%edge(ee,4) == -2) then 
                        volume_mesh%edge(ee,4) = cm2dopt%BC_zone_bc(ii)
                    end if 
                end if
            end if
        end do 
    end if 
end do 
return 
end subroutine set_custom_bcs




!Subroutine to remove far field adjacent regions ===========================
subroutine remove_farfield_adjacent(volume_mesh)
implicit none 

!Variables - Import
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ee,cc,vv,Nupdate,r3,r4
integer(in) :: Ncell_N,Nvtx_N,Nedge_N,Ncell_0,Nedge_0,Nvtx_0
integer(in) :: cell_remove(volume_mesh%ncell)
integer(in) :: edge_remove(volume_mesh%nedge)
integer(in) :: vtx_remove(volume_mesh%nvtx)
integer(in) :: cell_idx_N(volume_mesh%ncell)
integer(in) :: edge_idx_N(volume_mesh%nedge)
integer(in) :: vtx_idx_N(volume_mesh%nvtx)
integer(in) :: vtx_surfseg_temp(volume_mesh%nvtx)
integer(in) :: cell_level_temp(volume_mesh%ncell)
integer(in) :: edges_temp(volume_mesh%nedge,4)
real(dp) :: vtx_temp(volume_mesh%nvtx,2)

!Tag cells to remove as any with a far field boundary condition 
cell_remove(:) = 0 
do ee=1,volume_mesh%nedge
    if (volume_mesh%edge(ee,3) == -2) then 
        cell_remove(volume_mesh%edge(ee,4)) = 1
    elseif (volume_mesh%edge(ee,4) == -2) then 
        cell_remove(volume_mesh%edge(ee,3)) = 1
    end if
end do 

!Flood through mesh 
do cc=1,volume_mesh%ncell
    Nupdate = 0 
    do ee=1,volume_mesh%nedge
        if ((volume_mesh%edge(ee,3) .GT. 0) .AND. (volume_mesh%edge(ee,4) .GT. 0)) then 
            r3 = cell_remove(volume_mesh%edge(ee,3))
            r4= cell_remove(volume_mesh%edge(ee,4))
            if ((r3 == 1) .AND. (r4 == 0)) then
                cell_remove(volume_mesh%edge(ee,4)) = 1
                Nupdate = Nupdate + 1
            elseif ((r3 == 0) .AND. (r4 == 1)) then
                cell_remove(volume_mesh%edge(ee,3)) = 1
                Nupdate = Nupdate + 1
            end if 
        end if
    end do 
    if (Nupdate == 0) then 
        exit 
    end if 
end do 

!Tag edges to remove 
edge_remove(:) = 0 
do ee=1,volume_mesh%nedge
    if (volume_mesh%edge(ee,3) == -2) then 
        edge_remove(ee) = 1
    elseif (volume_mesh%edge(ee,4) == -2) then 
        edge_remove(ee) = 1
    elseif (volume_mesh%edge(ee,3) .GT. 0) then 
        if (cell_remove(volume_mesh%edge(ee,3)) == 1) then 
            edge_remove(ee) = 1
        end if 
    elseif (volume_mesh%edge(ee,4) .GT. 0) then 
        if (cell_remove(volume_mesh%edge(ee,4)) == 1) then 
            edge_remove(ee) = 1
        end if 
    end if
end do 

!Tag vertices to remove 
vtx_remove(:) = 0 
do ee=1,volume_mesh%nedge
    if (edge_remove(ee) == 1) then 
        vtx_remove(volume_mesh%edge(ee,1:2)) = 1
    end if 
end do 

!Re-index retained mesh region 
Ncell_N = 0 
cell_idx_N(:) = 0 
do cc=1,volume_mesh%ncell
    if (cell_remove(cc) == 0) then 
        Ncell_N = Ncell_N + 1
        cell_idx_N(cc) = Ncell_N
    end if
end do 
Nedge_N = 0 
edge_idx_N(:) = 0 
do ee=1,volume_mesh%nedge
    if (edge_remove(ee) == 0) then 
        Nedge_N = Nedge_N + 1
        edge_idx_N(ee) = Nedge_N
    end if
end do 
Nvtx_N = 0 
vtx_idx_N(:) = 0 
do vv=1,volume_mesh%nvtx
    if (vtx_remove(vv) == 0) then 
        Nvtx_N = Nvtx_N + 1
        vtx_idx_N(vv) = Nvtx_N
    end if
end do 

!Build new mesh 
Ncell_0 = volume_mesh%ncell
Nedge_0 = volume_mesh%nedge 
Nvtx_0 = volume_mesh%nvtx
cell_level_temp(:) = volume_mesh%cell_level(:)
edges_temp(:,:) = volume_mesh%edge(:,:)
vtx_temp(:,:) = volume_mesh%vertices(:,:)
vtx_surfseg_temp(:) = volume_mesh%vtx_surfseg(:)
deallocate(volume_mesh%cell_level)
deallocate(volume_mesh%edge)
deallocate(volume_mesh%vertices)
deallocate(volume_mesh%vtx_surfseg)
volume_mesh%ncell = Ncell_N
volume_mesh%nedge = Nedge_N
volume_mesh%nvtx = Nvtx_N
allocate(volume_mesh%edge(volume_mesh%nedge,4)) 
allocate(volume_mesh%vertices(volume_mesh%nvtx,2))
allocate(volume_mesh%vtx_surfseg(volume_mesh%nvtx))
allocate(volume_mesh%cell_level(volume_mesh%ncell))
do ee=1,Nedge_0
    if (edge_idx_N(ee) .NE. 0) then
        volume_mesh%edge(edge_idx_N(ee),:) = edges_temp(ee,:)
        volume_mesh%edge(edge_idx_N(ee),1) = vtx_idx_N(volume_mesh%edge(edge_idx_N(ee),1))
        volume_mesh%edge(edge_idx_N(ee),2) = vtx_idx_N(volume_mesh%edge(edge_idx_N(ee),2))
        if (volume_mesh%edge(edge_idx_N(ee),3) .GT. 0) then 
            volume_mesh%edge(edge_idx_N(ee),3) = cell_idx_N(volume_mesh%edge(edge_idx_N(ee),3))
        end if
        if (volume_mesh%edge(edge_idx_N(ee),4) .GT. 0) then 
            volume_mesh%edge(edge_idx_N(ee),4) = cell_idx_N(volume_mesh%edge(edge_idx_N(ee),4))
        end if
    end if 
end do 
do vv=1,Nvtx_0
    if (vtx_idx_N(vv) .NE. 0) then
        volume_mesh%vertices(vtx_idx_N(vv),:) = vtx_temp(vv,:)
        volume_mesh%vtx_surfseg(vtx_idx_N(vv)) = vtx_surfseg_temp(vv)
    end if 
end do  
do cc=1,Ncell_0
    if (cell_idx_N(cc) .NE. 0) then
        volume_mesh%cell_level(cell_idx_N(cc)) = cell_level_temp(cc)
    end if 
end do 
return 
end subroutine remove_farfield_adjacent




!Subroutine to remove non custom boundary condition adjacent regions ===========================
subroutine remove_ncb_adjacent(volume_mesh,cm2dopt)
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ee,ii,cc,vv,Nupdate,r3,r4,is_in_custom_zone
integer(in) :: Ncell_N,Nvtx_N,Nedge_N,Ncell_0,Nedge_0,Nvtx_0
integer(in) :: cell_remove(volume_mesh%ncell)
integer(in) :: edge_remove(volume_mesh%nedge)
integer(in) :: vtx_remove(volume_mesh%nvtx)
integer(in) :: cell_idx_N(volume_mesh%ncell)
integer(in) :: edge_idx_N(volume_mesh%nedge)
integer(in) :: vtx_idx_N(volume_mesh%nvtx)
integer(in) :: vtx_surfseg_temp(volume_mesh%nvtx)
integer(in) :: cell_level_temp(volume_mesh%ncell)
integer(in) :: edges_temp(volume_mesh%nedge,4)
real(dp) :: xmin,xmax,ymin,ymax
real(dp) :: vtx_temp(volume_mesh%nvtx,2)

!Tag cells to remove as any with a not custom set boundary condition 
cell_remove(:) = 0 
do ee=1,volume_mesh%nedge
    if ((volume_mesh%edge(ee,3) .LT. -1) .OR. (volume_mesh%edge(ee,4) .LT. -1)) then 

        !Edge coordinate bounds
        xmin = minval(volume_mesh%vertices(volume_mesh%edge(ee,1:2),1))
        xmax = maxval(volume_mesh%vertices(volume_mesh%edge(ee,1:2),1))
        ymin = minval(volume_mesh%vertices(volume_mesh%edge(ee,1:2),2))
        ymax = maxval(volume_mesh%vertices(volume_mesh%edge(ee,1:2),2))

        !Check if edge is within a custom boundary condition zone 
        is_in_custom_zone = 0
        do ii=1,cm2dopt%Nzone_cBC
            if ((xmin .LE. cm2dopt%BC_zone_coords(ii,2)) .AND. (xmax .GE. cm2dopt%BC_zone_coords(ii,1))) then 
                if ((ymin .LE. cm2dopt%BC_zone_coords(ii,4)) .AND. (ymax .GE. cm2dopt%BC_zone_coords(ii,3))) then 
                    is_in_custom_zone = 1
                    exit
                end if
            end if
        end do 

        !If not within a custom zone then tag adjacent cells to remove 
        if (is_in_custom_zone == 0) then 
            if (volume_mesh%edge(ee,3) .LT. 0) then 
                cell_remove(volume_mesh%edge(ee,4)) = 1
            elseif (volume_mesh%edge(ee,4) .LT. 0) then 
                cell_remove(volume_mesh%edge(ee,3)) = 1
            end if
        end if 
    end if
end do 

!Debug write remove tagged cells at the first itteration
! open(11,file='io/cell_remove')
! do cc=1,volume_mesh%ncell
!     write(11,*) cell_remove(cc)
! end do 
! close(11)

!Flood through mesh 
do cc=1,volume_mesh%ncell
    Nupdate = 0 
    do ee=1,volume_mesh%nedge
        if ((volume_mesh%edge(ee,3) .GT. 0) .AND. (volume_mesh%edge(ee,4) .GT. 0)) then 
            r3 = cell_remove(volume_mesh%edge(ee,3))
            r4= cell_remove(volume_mesh%edge(ee,4))
            if ((r3 == 1) .AND. (r4 == 0)) then
                cell_remove(volume_mesh%edge(ee,4)) = 1
                Nupdate = Nupdate + 1
            elseif ((r3 == 0) .AND. (r4 == 1)) then
                cell_remove(volume_mesh%edge(ee,3)) = 1
                Nupdate = Nupdate + 1
            end if 
        end if
    end do 
    if (Nupdate == 0) then 
        exit 
    end if 
end do 

!Tag edges to remove that are a part of a cell being removed 
edge_remove(:) = 0 
do ee=1,volume_mesh%nedge
    if (volume_mesh%edge(ee,3) .GT. 0) then 
        if (cell_remove(volume_mesh%edge(ee,3)) == 1) then 
            edge_remove(ee) = 1
        end if 
    end if 
    if (volume_mesh%edge(ee,4) .GT. 0) then 
        if (cell_remove(volume_mesh%edge(ee,4)) == 1) then 
            edge_remove(ee) = 1
        end if 
    end if 
end do 

!Tag vertices to remove 
vtx_remove(:) = 0 
do ee=1,volume_mesh%nedge
    if (edge_remove(ee) == 1) then 
        vtx_remove(volume_mesh%edge(ee,1:2)) = 1
    end if 
end do 

!Re-index retained mesh region 
Ncell_N = 0 
cell_idx_N(:) = 0 
do cc=1,volume_mesh%ncell
    if (cell_remove(cc) == 0) then 
        Ncell_N = Ncell_N + 1
        cell_idx_N(cc) = Ncell_N
    end if
end do 
Nedge_N = 0 
edge_idx_N(:) = 0 
do ee=1,volume_mesh%nedge
    if (edge_remove(ee) == 0) then 
        Nedge_N = Nedge_N + 1
        edge_idx_N(ee) = Nedge_N
    end if
end do 
Nvtx_N = 0 
vtx_idx_N(:) = 0 
do vv=1,volume_mesh%nvtx
    if (vtx_remove(vv) == 0) then 
        Nvtx_N = Nvtx_N + 1
        vtx_idx_N(vv) = Nvtx_N
    end if
end do 

!Build new mesh 
Ncell_0 = volume_mesh%ncell
Nedge_0 = volume_mesh%nedge 
Nvtx_0 = volume_mesh%nvtx
cell_level_temp(:) = volume_mesh%cell_level(:)
edges_temp(:,:) = volume_mesh%edge(:,:)
vtx_temp(:,:) = volume_mesh%vertices(:,:)
vtx_surfseg_temp(:) = volume_mesh%vtx_surfseg(:)
deallocate(volume_mesh%cell_level)
deallocate(volume_mesh%edge)
deallocate(volume_mesh%vertices)
deallocate(volume_mesh%vtx_surfseg)
volume_mesh%ncell = Ncell_N
volume_mesh%nedge = Nedge_N
volume_mesh%nvtx = Nvtx_N
allocate(volume_mesh%edge(volume_mesh%nedge,4)) 
allocate(volume_mesh%vertices(volume_mesh%nvtx,2))
allocate(volume_mesh%vtx_surfseg(volume_mesh%nvtx))
allocate(volume_mesh%cell_level(volume_mesh%ncell))
do ee=1,Nedge_0
    if (edge_idx_N(ee) .NE. 0) then
        volume_mesh%edge(edge_idx_N(ee),:) = edges_temp(ee,:)
        volume_mesh%edge(edge_idx_N(ee),1) = vtx_idx_N(volume_mesh%edge(edge_idx_N(ee),1))
        volume_mesh%edge(edge_idx_N(ee),2) = vtx_idx_N(volume_mesh%edge(edge_idx_N(ee),2))
        if (volume_mesh%edge(edge_idx_N(ee),3) .GT. 0) then 
            volume_mesh%edge(edge_idx_N(ee),3) = cell_idx_N(volume_mesh%edge(edge_idx_N(ee),3))
        end if
        if (volume_mesh%edge(edge_idx_N(ee),4) .GT. 0) then 
            volume_mesh%edge(edge_idx_N(ee),4) = cell_idx_N(volume_mesh%edge(edge_idx_N(ee),4))
        end if
    end if 
end do 
do vv=1,Nvtx_0
    if (vtx_idx_N(vv) .NE. 0) then
        volume_mesh%vertices(vtx_idx_N(vv),:) = vtx_temp(vv,:)
        volume_mesh%vtx_surfseg(vtx_idx_N(vv)) = vtx_surfseg_temp(vv)
    end if 
end do  
do cc=1,Ncell_0
    if (cell_idx_N(cc) .NE. 0) then
        volume_mesh%cell_level(cell_idx_N(cc)) = cell_level_temp(cc)
    end if 
end do 
return 
end subroutine remove_ncb_adjacent




!Subroutine to remove isolated wall adjacent only regions ===========================
subroutine remove_isolated_regions(volume_mesh)
implicit none 

!Variables - Import
type(vol_mesh_data) :: volume_mesh


!Variables - Local 
integer(in) :: ee,cc,vv,Nupdate,r3,r4
integer(in) :: Ncell_N,Nvtx_N,Nedge_N,Ncell_0,Nedge_0,Nvtx_0
integer(in) :: cell_keep(volume_mesh%ncell)
integer(in) :: cell_remove(volume_mesh%ncell)
integer(in) :: edge_remove(volume_mesh%nedge)
integer(in) :: vtx_remove(volume_mesh%nvtx)
integer(in) :: cell_idx_N(volume_mesh%ncell)
integer(in) :: edge_idx_N(volume_mesh%nedge)
integer(in) :: vtx_idx_N(volume_mesh%nvtx)
integer(in) :: vtx_surfseg_temp(volume_mesh%nvtx)
integer(in) :: cell_level_temp(volume_mesh%ncell)
integer(in) :: edges_temp(volume_mesh%nedge,4)
real(dp) :: vtx_temp(volume_mesh%nvtx,2)


!Tag initial cells to keep as any with any form of entry/exit/farfield boundary condition 
cell_keep(:) = 0 
do ee=1,volume_mesh%nedge
    if (volume_mesh%edge(ee,3) .LE. -2) then 
        cell_keep(volume_mesh%edge(ee,4)) = 1
    elseif (volume_mesh%edge(ee,4) .LE. -2) then 
        cell_keep(volume_mesh%edge(ee,3)) = 1
    end if
end do 

!Flood through mesh 
do cc=1,volume_mesh%ncell
    Nupdate = 0 
    do ee=1,volume_mesh%nedge
        if ((volume_mesh%edge(ee,3) .GT. 0) .AND. (volume_mesh%edge(ee,4) .GT. 0)) then 
            r3 = cell_keep(volume_mesh%edge(ee,3))
            r4 = cell_keep(volume_mesh%edge(ee,4))
            if ((r3 == 1) .AND. (r4 == 0)) then
                cell_keep(volume_mesh%edge(ee,4)) = 1
                Nupdate = Nupdate + 1
            elseif ((r3 == 0) .AND. (r4 == 1)) then
                cell_keep(volume_mesh%edge(ee,3)) = 1
                Nupdate = Nupdate + 1
            end if 
        end if
    end do 
    if (Nupdate == 0) then 
        exit 
    end if 
end do 

!Tag cells to remove as any not tagged to keep 
cell_remove(:) = 0 
do cc=1,volume_mesh%ncell
    if (cell_keep(cc) == 0) then 
        cell_remove(cc) = 1 
    end if 
end do 

!Tag edges to remove 
edge_remove(:) = 0 
do ee=1,volume_mesh%nedge
    if (volume_mesh%edge(ee,3) .GT. 0) then 
        if (cell_remove(volume_mesh%edge(ee,3)) == 1) then 
            edge_remove(ee) = 1
        end if 
    elseif (volume_mesh%edge(ee,4) .GT. 0) then 
        if (cell_remove(volume_mesh%edge(ee,4)) == 1) then 
            edge_remove(ee) = 1
        end if 
    end if
end do 

!Tag vertices to remove 
vtx_remove(:) = 0 
do ee=1,volume_mesh%nedge
    if (edge_remove(ee) == 1) then 
        vtx_remove(volume_mesh%edge(ee,1:2)) = 1
    end if 
end do 

!Re-index retained mesh region 
Ncell_N = 0 
cell_idx_N(:) = 0 
do cc=1,volume_mesh%ncell
    if (cell_remove(cc) == 0) then 
        Ncell_N = Ncell_N + 1
        cell_idx_N(cc) = Ncell_N
    end if
end do 
Nedge_N = 0 
edge_idx_N(:) = 0 
do ee=1,volume_mesh%nedge
    if (edge_remove(ee) == 0) then 
        Nedge_N = Nedge_N + 1
        edge_idx_N(ee) = Nedge_N
    end if
end do 
Nvtx_N = 0 
vtx_idx_N(:) = 0 
do vv=1,volume_mesh%nvtx
    if (vtx_remove(vv) == 0) then 
        Nvtx_N = Nvtx_N + 1
        vtx_idx_N(vv) = Nvtx_N
    end if
end do 

!Build new mesh 
Ncell_0 = volume_mesh%ncell
Nedge_0 = volume_mesh%nedge 
Nvtx_0 = volume_mesh%nvtx
cell_level_temp(:) = volume_mesh%cell_level(:)
edges_temp(:,:) = volume_mesh%edge(:,:)
vtx_temp(:,:) = volume_mesh%vertices(:,:)
vtx_surfseg_temp(:) = volume_mesh%vtx_surfseg(:)
deallocate(volume_mesh%cell_level)
deallocate(volume_mesh%edge)
deallocate(volume_mesh%vertices)
deallocate(volume_mesh%vtx_surfseg)
volume_mesh%ncell = Ncell_N
volume_mesh%nedge = Nedge_N
volume_mesh%nvtx = Nvtx_N
allocate(volume_mesh%edge(volume_mesh%nedge,4)) 
allocate(volume_mesh%vertices(volume_mesh%nvtx,2))
allocate(volume_mesh%vtx_surfseg(volume_mesh%nvtx))
allocate(volume_mesh%cell_level(volume_mesh%ncell))
do ee=1,Nedge_0
    if (edge_idx_N(ee) .NE. 0) then
        volume_mesh%edge(edge_idx_N(ee),:) = edges_temp(ee,:)
        volume_mesh%edge(edge_idx_N(ee),1) = vtx_idx_N(volume_mesh%edge(edge_idx_N(ee),1))
        volume_mesh%edge(edge_idx_N(ee),2) = vtx_idx_N(volume_mesh%edge(edge_idx_N(ee),2))
        if (volume_mesh%edge(edge_idx_N(ee),3) .GT. 0) then 
            volume_mesh%edge(edge_idx_N(ee),3) = cell_idx_N(volume_mesh%edge(edge_idx_N(ee),3))
        end if
        if (volume_mesh%edge(edge_idx_N(ee),4) .GT. 0) then 
            volume_mesh%edge(edge_idx_N(ee),4) = cell_idx_N(volume_mesh%edge(edge_idx_N(ee),4))
        end if
    end if 
end do 
do vv=1,Nvtx_0
    if (vtx_idx_N(vv) .NE. 0) then
        volume_mesh%vertices(vtx_idx_N(vv),:) = vtx_temp(vv,:)
        volume_mesh%vtx_surfseg(vtx_idx_N(vv)) = vtx_surfseg_temp(vv)
    end if 
end do  
do cc=1,Ncell_0
    if (cell_idx_N(cc) .NE. 0) then
        volume_mesh%cell_level(cell_idx_N(cc)) = cell_level_temp(cc)
    end if 
end do 
return 
end subroutine remove_isolated_regions




!Subroutine to flip edges with a target boundary condition ===========================
subroutine flip_set_edges(volume_mesh,fcondition)
implicit none 

!Variables - Import
integer(in) :: fcondition
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ee,c3,c4,v1,v2

!Flip edges with the target condition fcondition as cl or cr 
do ee=1,volume_mesh%nedge
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    if ((c3 == fcondition) .OR. (c4 == fcondition)) then 
        v1 = volume_mesh%edge(ee,1)
        v2 = volume_mesh%edge(ee,2)
        volume_mesh%edge(ee,1) = v2
        volume_mesh%edge(ee,2) = v1
        volume_mesh%edge(ee,3) = c4
        volume_mesh%edge(ee,4) = c3
    end if 
end do 
return 
end subroutine flip_set_edges




!Subroutine to remove double boundary condition edges ===========================
subroutine remove_double_bc_edges(volume_mesh,cm2dopt)
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ee,vv,nremove,NvtxN,NedgeN,eins 
integer(in) :: eremove(volume_mesh%nedge),vremove(volume_mesh%nvtx)
integer(in) :: vtx_map(volume_mesh%nvtx),vtx_surfseg_temp(volume_mesh%nvtx),etemp(volume_mesh%nedge,4) 
real(dp) :: vtemp(volume_mesh%nvtx,2)

!Tag edges to remove 
nremove = 0 
eremove(:) = 0 
do ee=1,volume_mesh%nedge 
    if ((volume_mesh%edge(ee,3) .LE. 0) .AND. (volume_mesh%edge(ee,4) .LE. 0)) then 
        nremove = nremove + 1
        eremove(ee) = 1
    end if 
end do 

!Tag vertices to remove 
vremove(:) = 0 
do ee=1,volume_mesh%nedge 
    if (eremove(ee) == 1) then 
        vremove(volume_mesh%edge(ee,1:2)) = 1
    end if 
end do 
do ee=1,volume_mesh%nedge 
    if (eremove(ee) == 0) then 
        vremove(volume_mesh%edge(ee,1:2)) = 0
    end if 
end do 

!Map vertices and edges 
NvtxN = 0 
vtx_map(:) = 0 
do vv=1,volume_mesh%nvtx
    if (vremove(vv) == 0) then
        NvtxN = NvtxN + 1
        vtx_map(vv) = NvtxN
    end if 
end do 
NedgeN = 0 
do ee=1,volume_mesh%nedge 
    if (eremove(ee) == 0) then 
        NedgeN = NedgeN + 1
    end if
end do 

!Reconstruct mesh arrays 
etemp(:,:) = volume_mesh%edge(:,:)
vtemp(:,:) = volume_mesh%vertices(:,:)
vtx_surfseg_temp(:) = volume_mesh%vtx_surfseg(:)
deallocate(volume_mesh%edge)
deallocate(volume_mesh%vertices)
deallocate(volume_mesh%vtx_surfseg)
allocate(volume_mesh%vertices(NvtxN,2))
allocate(volume_mesh%edge(NedgeN,4))
allocate(volume_mesh%vtx_surfseg(NvtxN))
eins = 0 
do ee=1,volume_mesh%nedge 
    if (eremove(ee) == 0) then
        eins = eins + 1
        volume_mesh%edge(eins,1) = vtx_map(etemp(ee,1))
        volume_mesh%edge(eins,2) = vtx_map(etemp(ee,2))
        volume_mesh%edge(eins,3:4) = etemp(ee,3:4)
    end if 
end do 
do vv=1,volume_mesh%nvtx 
    if (vremove(vv) == 0) then
        volume_mesh%vertices(vtx_map(vv),:) = vtemp(vv,:)
        volume_mesh%vtx_surfseg(vtx_map(vv)) = vtx_surfseg_temp(vv)
    end if
end do 
volume_mesh%nvtx = NvtxN
volume_mesh%nedge = NedgeN

!Display
if (nremove .GT. 0) then 
    if (cm2dopt%dispt == 1) then
        write(*,'(A,I0,A)') '    {removed ',nremove,' double boundary condition edges}'
    end if 
end if 
return 
end subroutine remove_double_bc_edges




!Subroutine to simplify the mesh surface ===========================
subroutine simplify_surface(volume_mesh)
implicit none 

!Variables - Import
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ee,vv,cc,ff,aa
integer(in) :: nedgeN,nvtxN,v1,v2,ebase,etgt,nupdate,eadj,emtgt,vstart,vend,nvem,vtgt,vadj,vp,vc,eposfound
integer(in) :: edge_map(volume_mesh%nedge),vtx_map(volume_mesh%nvtx)
integer(in) :: edge_merge_idx(volume_mesh%nedge),cell_nsedge(volume_mesh%ncell),medge_evtxs(volume_mesh%nvtx)
integer(in) :: edgeN(volume_mesh%nedge,4),v2e_vsurf(volume_mesh%nvtx,2)
integer(in), dimension(:), allocatable ::  vtx_surfsegN
integer(in), dimension(:,:), allocatable :: cell_surfedges
real(dp) :: distc,distp,ef,emx,emy
real(dp) :: medge_evtxs_dfrac(volume_mesh%nvtx)
real(dp) :: edgemidpN(volume_mesh%nedge,2)
real(dp), dimension(:,:), allocatable :: verticesN

!Index retained non surface edges 
nedgeN = 0 
edge_map(:) = 0 
do ee=1,volume_mesh%nedge
    if ((volume_mesh%edge(ee,3) .NE. -1) .AND. (volume_mesh%edge(ee,4) .NE. -1)) then
        nedgeN = nedgeN + 1
        edge_map(ee) = nedgeN
    end if 
end do 

!Copy the retained edges to the new strucrure
edgeN(:,:) = 0 
do ee=1,volume_mesh%nedge
    if (edge_map(ee) .NE. 0) then 
        edgeN(edge_map(ee),:) = volume_mesh%edge(ee,:)
    end if 
end do 

!Build v2e for surface vertices and surface edges 
v2e_vsurf(:,:) = 0 
do ee=1,volume_mesh%nedge
    if ((volume_mesh%edge(ee,3) == -1) .OR. (volume_mesh%edge(ee,4) == -1)) then
        v2e_vsurf(volume_mesh%edge(ee,1),2) = ee 
        v2e_vsurf(volume_mesh%edge(ee,2),1) = ee 
    end if
end do 

!Accumulate surface edges to each cell 
cell_nsedge(:) = 0 
do ee=1,volume_mesh%nedge
    if (volume_mesh%edge(ee,3) == -1) then
        cell_nsedge(volume_mesh%edge(ee,4)) = cell_nsedge(volume_mesh%edge(ee,4)) + 1
    elseif (volume_mesh%edge(ee,4) == -1) then 
        cell_nsedge(volume_mesh%edge(ee,3)) = cell_nsedge(volume_mesh%edge(ee,3)) + 1
    end if
end do 
allocate(cell_surfedges(volume_mesh%ncell,maxval(cell_nsedge(:))))
cell_surfedges(:,:) = 0 
cell_nsedge(:) = 0 
do ee=1,volume_mesh%nedge
    if (volume_mesh%edge(ee,3) == -1) then
        cell_nsedge(volume_mesh%edge(ee,4)) = cell_nsedge(volume_mesh%edge(ee,4)) + 1
        cell_surfedges(volume_mesh%edge(ee,4),cell_nsedge(volume_mesh%edge(ee,4))) = ee 
    elseif (volume_mesh%edge(ee,4) == -1) then 
        cell_nsedge(volume_mesh%edge(ee,3)) = cell_nsedge(volume_mesh%edge(ee,3)) + 1
        cell_surfedges(volume_mesh%edge(ee,3),cell_nsedge(volume_mesh%edge(ee,3))) = ee 
    end if
end do 

!Build simplified edges in surface cells 
edge_merge_idx(:) = 0 
do cc=1,volume_mesh%ncell
    if (cell_nsedge(cc) .GT. 0) then 
        do ff=1,volume_mesh%nedge

            !Find an edge in this cell that is not tagged with a new merge index
            ebase = 0 
            do ee=1,cell_nsedge(cc)
                etgt = cell_surfedges(cc,ee)
                if (edge_merge_idx(etgt) == 0) then 
                    ebase = etgt 
                    exit 
                end if 
            end do 

            !If all edges are tagged then exit as new edges are complete here 
            if (ebase == 0) then 
                exit 
            end if 

            !Increment edge count 
            nedgeN = nedgeN + 1

            !Set the tag of this selected edge ebase
            edge_merge_idx(ebase) = nedgeN

            !Set adjacency of this new edge nedgeN
            edgeN(nedgeN,3:4) = volume_mesh%edge(ebase,3:4)

            !Flood the new edge tag untill a cell boundary is reached 
            do aa=1,cell_nsedge(cc)
                nupdate = 0 
                do ee=1,cell_nsedge(cc)

                    !Edge
                    etgt = cell_surfedges(cc,ee)
                    v1 = volume_mesh%edge(etgt,1)
                    v2 = volume_mesh%edge(etgt,2)

                    !If non zero tagged 
                    if (edge_merge_idx(etgt) .NE. 0) then 
                    
                        !Flood tag across end 1 
                        eadj = 0 
                        if (v2e_vsurf(v1,1) == etgt) then 
                            eadj = v2e_vsurf(v1,2)
                        elseif (v2e_vsurf(v1,2) == etgt) then 
                            eadj = v2e_vsurf(v1,1)
                        end if 
                        if (eadj .NE. 0) then 
                            if ((volume_mesh%edge(eadj,3) == cc) .OR. (volume_mesh%edge(eadj,4) == cc)) then 
                                if (edge_merge_idx(eadj) == 0) then 
                                    edge_merge_idx(eadj) = edge_merge_idx(etgt)
                                    nupdate = nupdate + 1
                                end if 
                            end if 
                        end if 

                        !Flood tag across end 2 
                        eadj = 0 
                        if (v2e_vsurf(v2,1) == etgt) then 
                            eadj = v2e_vsurf(v2,2)
                        elseif (v2e_vsurf(v2,2) == etgt) then 
                            eadj = v2e_vsurf(v2,1)
                        end if 
                        if (eadj .NE. 0) then 
                            if ((volume_mesh%edge(eadj,3) == cc) .OR. (volume_mesh%edge(eadj,4) == cc)) then 
                                if (edge_merge_idx(eadj) == 0) then 
                                    edge_merge_idx(eadj) = edge_merge_idx(etgt)
                                    nupdate = nupdate + 1
                                end if 
                            end if 
                        end if 
                    end if 
                end do 
                if (nupdate == 0) then 
                    exit 
                end if 
            end do 
        end do 
    end if 
end do 

!Build new simplified edges 
do cc=1,volume_mesh%ncell
    if (cell_nsedge(cc) .GT. 0) then
        do ee=1,cell_nsedge(cc)

            !Existing edge
            etgt = cell_surfedges(cc,ee)

            !Edge ends
            v1 = volume_mesh%edge(etgt,1)
            v2 = volume_mesh%edge(etgt,2)

            !If this existing edge is an end of edge edge_merge_idx(etgt) add the required vertices to build new edges ---
            !End 1 
            eadj = 0 
            if (v2e_vsurf(v1,1) == etgt) then 
                eadj = v2e_vsurf(v1,2)
            elseif (v2e_vsurf(v1,2) == etgt) then 
                eadj = v2e_vsurf(v1,1)
            end if 
            if (eadj .NE. 0) then 
                if (edge_merge_idx(eadj) .NE. edge_merge_idx(etgt)) then !add vertex v1 as start of etgt and end of eadj
                    edgeN(edge_merge_idx(etgt),1) = v1
                    edgeN(edge_merge_idx(eadj),2) = v1
                end if 
            else !add vertex v1 as start of etgt
                edgeN(edge_merge_idx(etgt),1) = v1
            end if 

            !End 2 
            eadj = 0 
            if (v2e_vsurf(v2,1) == etgt) then 
                eadj = v2e_vsurf(v2,2)
            elseif (v2e_vsurf(v2,2) == etgt) then 
                eadj = v2e_vsurf(v2,1)
            end if 
            if (eadj .NE. 0) then 
                if (edge_merge_idx(eadj) .NE. edge_merge_idx(etgt)) then !add vertex v2 as end of etgt and start of eadj
                    edgeN(edge_merge_idx(etgt),2) = v2
                    edgeN(edge_merge_idx(eadj),1) = v2
                end if
            else !add vertex v2 as end of etgt
                edgeN(edge_merge_idx(etgt),2) = v2
            end if

            !Debug --
            !print *, v2e_vsurf(v1,:),' || ',v2e_vsurf(v2,:),' vtx: ',v1,' / ',v2
        end do 
    end if
end do 

!Find the on surface midpoint of each simplified edge 
medge_evtxs(:) = 0 
medge_evtxs_dfrac(:) = 0.0d0 
edgemidpN(:,:) = 0.0d0 
do cc=1,volume_mesh%ncell
    if (cell_nsedge(cc) .GT. 0) then
        do aa=1,2*cell_nsedge(cc)
            
            !Find merged edge 
            emtgt = 0 
            do ee=1,cell_nsedge(cc)
                etgt = cell_surfedges(cc,ee)
                if (edge_merge_idx(etgt) .NE. 0) then 
                    emtgt = edge_merge_idx(etgt)
                end if 
            end do 

            !Exit with no merged edge 
            if (emtgt == 0) then 
                exit 
            end if 

            !Order existing vertices along this edge 
            vstart = edgeN(emtgt,1)
            vend = edgeN(emtgt,1)
            nvem = 1
            medge_evtxs(1) = vstart
            vtgt = vstart
            do ee=1,2*cell_nsedge(cc)

                !Find next edge on the surface 
                eadj = 0 
                etgt = v2e_vsurf(vtgt,2)
                if (etgt .GT. 0) then 
                    if (edge_merge_idx(etgt) == emtgt) then 
                        eadj = etgt
                    end if 
                end if 

                !Exit at complete search 
                if (eadj == 0) then 
                    exit 
                end if 

                !Find next vertex
                if (volume_mesh%edge(eadj,1) == vtgt) then 
                    vadj = volume_mesh%edge(eadj,2)
                else
                    vadj = volume_mesh%edge(eadj,1)
                end if 

                !Store next vertex
                nvem = nvem + 1
                medge_evtxs(nvem) = vadj

                !Update base vertex
                vtgt = vadj

                !Exit if final vertex
                if (vadj == vend) then 
                    exit 
                end if 
            end do 

            !Find distance fraction at each vertex 
            do ee=2,nvem
                vp = medge_evtxs(ee-1)
                vc = medge_evtxs(ee)
                medge_evtxs_dfrac(ee) = medge_evtxs_dfrac(ee-1) + norm2(volume_mesh%vertices(vc,:) - volume_mesh%vertices(vp,:))
            end do 
            medge_evtxs_dfrac(1:nvem) = medge_evtxs_dfrac(1:nvem)/medge_evtxs_dfrac(nvem)

            !Find position at 0.5 distance fraction 
            eposfound = 0 
            do ee=2,nvem
                vp = ee-1
                vc = ee
                if ((medge_evtxs_dfrac(vc) .GE. 0.5d0) .AND. (medge_evtxs_dfrac(vp) .LE. 0.5d0)) then 
                    distc = medge_evtxs_dfrac(vc)
                    distp = medge_evtxs_dfrac(vp)
                    if ((distc - distp) == 0.0d0) then 
                        ef = 0.0d0 
                    else
                        ef = (0.5d0 - distp)/(distc - distp)
                    end if 
                    vp = medge_evtxs(ee-1)
                    vc = medge_evtxs(ee)
                    emx = volume_mesh%vertices(vc,1)*ef + volume_mesh%vertices(vp,1)*(1.0d0 - ef)
                    emy = volume_mesh%vertices(vc,2)*ef + volume_mesh%vertices(vp,2)*(1.0d0 - ef)
                    eposfound = 1
                    exit 
                end if 
            end do 
            if (eposfound == 1) then 
                edgemidpN(emtgt,1) = emx
                edgemidpN(emtgt,2) = emy
            else
                edgemidpN(emtgt,1) = 0.5d0*(volume_mesh%vertices(vstart,1) + volume_mesh%vertices(vend,1))
                edgemidpN(emtgt,2) = 0.5d0*(volume_mesh%vertices(vstart,2) + volume_mesh%vertices(vend,2))
            end if 

            !Reset vertex list and distance fractions 
            medge_evtxs(1:nvem) = 0
            medge_evtxs_dfrac(1:nvem) = 0.0d0 
        end do 
    end if
end do 

!Store new mesh edges and midpoints
volume_mesh%nedge = nedgeN !+ esnidx
deallocate(volume_mesh%edge)
allocate(volume_mesh%edge(volume_mesh%nedge,4))
volume_mesh%edge(:,:) = edgeN(1:nedgeN,:)
allocate(volume_mesh%edge_midpoint(volume_mesh%nedge,2))
volume_mesh%edge_midpoint(:,:) = edgemidpN(1:nedgeN,:)

!Map retained vertices 
nvtxN = 0 
vtx_map(:) = 0 
do ee=1,volume_mesh%nedge
    if (vtx_map(volume_mesh%edge(ee,1)) == 0) then 
        nvtxN = nvtxN + 1
        vtx_map(volume_mesh%edge(ee,1)) = nvtxN
    end if 
    if (vtx_map(volume_mesh%edge(ee,2)) == 0) then 
        nvtxN = nvtxN + 1
        vtx_map(volume_mesh%edge(ee,2)) = nvtxN
    end if 
end do 
do ee=1,volume_mesh%nedge
    volume_mesh%edge(ee,1) = vtx_map(volume_mesh%edge(ee,1))
    volume_mesh%edge(ee,2) = vtx_map(volume_mesh%edge(ee,2))
end do 

!Map retained vertices
allocate(verticesN(nvtxN,2))
do vv=1,volume_mesh%nvtx
    if (vtx_map(vv) .NE. 0) then 
        verticesN(vtx_map(vv),:) = volume_mesh%vertices(vv,:)
    end if 
end do 

!Update surface vertex links volume_mesh%vtx_surfseg
allocate(vtx_surfsegN(volume_mesh%nvtx))
do vv=1,volume_mesh%nvtx
    if (vtx_map(vv) .NE. 0) then 
        vtx_surfsegN(vtx_map(vv)) = volume_mesh%vtx_surfseg(vv)
    end if 
end do 
deallocate(volume_mesh%vtx_surfseg)
allocate(volume_mesh%vtx_surfseg(volume_mesh%nvtx))
volume_mesh%vtx_surfseg(:) = vtx_surfsegN(:)

!Store retained vertices
volume_mesh%nvtx = nvtxN
deallocate(volume_mesh%vertices)
allocate(volume_mesh%vertices(volume_mesh%nvtx,2))
volume_mesh%vertices(:,:) = verticesN(:,:)
return 
end subroutine simplify_surface




!Subroutine to remove mesh internal valence two vertices by merging adjacent edges ===========================
subroutine clean_internal_vlnc2_vertices(volume_mesh,cm2dopt)
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Local
integer(in) :: cc,ee,vv,aa 
integer(in) :: v1,v2,c3,c4,e1,e2,nintv2,Nenew,Nvnew,maxcedge,etgt,nmerge,emidx,vm1,vm2,c3n,c4n
integer(in) :: valence(volume_mesh%nvtx),vtx_bndry(volume_mesh%nvtx),edge_map(volume_mesh%nedge)
integer(in) :: nedge_cell(volume_mesh%ncell),vtx_map(volume_mesh%nvtx)
integer(in) :: V2E(volume_mesh%nvtx,4),V2Eloc(volume_mesh%nvtx,4)
integer(in), dimension(:), allocatable :: cedge_loop,cvtx_loop,emerge,vtx_surfseg_new
integer(in), dimension(:,:), allocatable :: cell2edge,edge_new
real(dp), dimension(:,:), allocatable :: vertices_new

!Evaluate vertex valence 
valence(:) = 0 
do ee=1,volume_mesh%nedge
    valence(volume_mesh%edge(ee,1)) = valence(volume_mesh%edge(ee,1)) + 1
    valence(volume_mesh%edge(ee,2)) = valence(volume_mesh%edge(ee,2)) + 1
end do 

!Tag surface vertices 
vtx_bndry(:) = 0
do ee=1,volume_mesh%nedge !all surface and boundary edges are ordered 
    v1 = volume_mesh%edge(ee,1)
    v2 = volume_mesh%edge(ee,2)
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    if ((c3 .LT. 0) .OR. (c4 .LT. 0)) then
        vtx_bndry(v1) = 1
        vtx_bndry(v2) = 1
    end if
end do 

!Count vertices that will be removed
nintv2 = 0 
do vv=1,volume_mesh%nvtx
    if ((vtx_bndry(vv) == 0) .AND. (valence(vv) == 2)) then 
        nintv2 = nintv2 + 1
    end if
end do 

!Build list of edges in each cell 
nedge_cell(:) = 0 
do ee=1,volume_mesh%nedge
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    if (c3 .GT. 0) then 
        nedge_cell(c3) = nedge_cell(c3) + 1
    end if 
    if (c4 .GT. 0) then 
        nedge_cell(c4) = nedge_cell(c4) + 1
    end if 
end do 
maxcedge = maxval(nedge_cell)
allocate(cell2edge(volume_mesh%ncell,maxcedge))
cell2edge(:,:) = 0 
nedge_cell(:) = 0 
do ee=1,volume_mesh%nedge
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    if (c3 .GT. 0) then 
        nedge_cell(c3) = nedge_cell(c3) + 1
        cell2edge(c3,nedge_cell(c3)) = ee 
    end if 
    if (c4 .GT. 0) then 
        nedge_cell(c4) = nedge_cell(c4) + 1
        cell2edge(c4,nedge_cell(c4)) = ee 
    end if 
end do 

!Build V2E on the mesh
V2E(:,:) = 0 
do ee=1,volume_mesh%nedge
    v1 = volume_mesh%edge(ee,1)
    v2 = volume_mesh%edge(ee,2)
    do aa=1,4
        if (V2E(v1,aa) == 0) then 
            V2E(v1,aa) = ee 
            exit 
        end if 
    end do 
    do aa=1,4
        if (V2E(v2,aa) == 0) then 
            V2E(v2,aa) = ee 
            exit 
        end if 
    end do 
end do 

!Map edges to merged indecies 
Nenew = 0 
edge_map(:) = 0 
do ee=1,volume_mesh%nedge
    v1 = volume_mesh%edge(ee,1)
    v2 = volume_mesh%edge(ee,2)
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    if ((c3 .LT. 0) .OR. (c4 .LT. 0)) then 
        if (edge_map(ee) == 0) then 
            Nenew = Nenew + 1
            edge_map(ee) = Nenew
        end if 
    elseif (valence(v1) == 2) then 
        e1 = V2E(v1,1)
        e2 = V2E(v1,2)
        if ((edge_map(e1) == 0) .AND. (edge_map(e2) == 0)) then 
            Nenew = Nenew + 1
            edge_map(e1) = -Nenew
            edge_map(e2) = -Nenew
        elseif ((edge_map(e1) == 0) .AND. (edge_map(e2) .NE. 0)) then 
            edge_map(e1) = edge_map(e2)
        elseif ((edge_map(e1) .NE. 0) .AND. (edge_map(e2) == 0)) then 
            edge_map(e2) = edge_map(e1)
        end if 
    elseif (valence(v2) == 2) then 
        e1 = V2E(v2,1)
        e2 = V2E(v2,2)
        if ((edge_map(e1) == 0) .AND. (edge_map(e2) == 0)) then 
            Nenew = Nenew + 1
            edge_map(e1) = -Nenew
            edge_map(e2) = -Nenew
        elseif ((edge_map(e1) == 0) .AND. (edge_map(e2) .NE. 0)) then 
            edge_map(e1) = edge_map(e2)
        elseif ((edge_map(e1) .NE. 0) .AND. (edge_map(e2) == 0)) then 
            edge_map(e2) = edge_map(e1)
        end if 
    else
        if (edge_map(ee) == 0) then 
            Nenew = Nenew + 1
            edge_map(ee) = Nenew
        end if 
    end if 
end do 

!Build new merged edges 
allocate(cedge_loop(maxcedge)) 
allocate(cvtx_loop(maxcedge))
allocate(emerge(maxcedge))
allocate(edge_new(Nenew,4))
edge_new(:,:) = 0 
V2Eloc(:,:) = 0 
do ee=1,volume_mesh%nedge
    if (edge_map(ee) .GT. 0) then 
        edge_new(edge_map(ee),:) = volume_mesh%edge(ee,:)
    end if 
end do 
do cc=1,volume_mesh%ncell

    !Build ordered loop in this cell 
    call build_cell_edges_vertices(cedge_loop,cvtx_loop,volume_mesh,nedge_cell,cell2edge,V2Eloc,cc)

    !Build merged edges
    do aa=1,nedge_cell(cc)

        !Find initial edge to merge 
        emidx = 0
        do ee=1,nedge_cell(cc)
            etgt = cell2edge(cc,ee)
            if (edge_map(etgt) .LT. 0) then 
                emidx = edge_map(etgt)
                exit 
            end if 
        end do 

        !Exit if no edges to merge 
        if (emidx == 0) then 
            exit 
        end if 

        !Collect edges  
        nmerge = 0 
        emerge(:) = 0 
        do ee=1,nedge_cell(cc)
            etgt = cell2edge(cc,ee)
            if (edge_map(etgt) == emidx) then 
                nmerge = nmerge + 1
                emerge(nmerge) = etgt 
                edge_map(etgt) = abs(edge_map(etgt))
            end if 
        end do 

        !Find start vertex and cell adjacency for this edge 
        etgt = emerge(1)
        if (valence(volume_mesh%edge(etgt,1)) .GT. 2) then 
            vm1 = volume_mesh%edge(etgt,1)
            c3n = volume_mesh%edge(etgt,3)
            c4n = volume_mesh%edge(etgt,4)
        elseif (valence(volume_mesh%edge(etgt,2)) .GT. 2) then 
            vm1 = volume_mesh%edge(etgt,2)
            c3n = volume_mesh%edge(etgt,4)
            c4n = volume_mesh%edge(etgt,3)
        else
            vm1 = 0 
            print *, '** edge merge failure'
        end if 

        !Find end vertex
        etgt = emerge(nmerge)
        if (valence(volume_mesh%edge(etgt,1)) .GT. 2) then 
            vm2 = volume_mesh%edge(etgt,1)
        elseif (valence(volume_mesh%edge(etgt,2)) .GT. 2) then 
            vm2 = volume_mesh%edge(etgt,2)
        else
            vm2 = 0 
            print *, '** edge merge failure'
        end if 

        !Build new edge 
        edge_new(abs(emidx),1) = vm1
        edge_new(abs(emidx),2) = vm2
        edge_new(abs(emidx),3) = c3n
        edge_new(abs(emidx),4) = c4n 
    end do 

    !Reset V2Eloc
    do ee=1,nedge_cell(cc)
        etgt = cell2edge(cc,ee)
        v1 = volume_mesh%edge(etgt,1)
        v2 = volume_mesh%edge(etgt,2)
        V2Eloc(v1,:) = 0 
        V2Eloc(v2,:) = 0
    end do 
end do 

!Build new vertices 
Nvnew = 0 
vtx_map(:) = 0 
do ee=1,Nenew
    v1 = edge_new(ee,1)
    v2 = edge_new(ee,2)
    if (vtx_map(v1) == 0) then 
        Nvnew = Nvnew + 1
        vtx_map(v1) = Nvnew
    end if 
    if (vtx_map(v2) == 0) then 
        Nvnew = Nvnew + 1
        vtx_map(v2) = Nvnew
    end if 
end do 
allocate(vertices_new(Nvnew,2))
allocate(vtx_surfseg_new(Nvnew))
vertices_new(:,:) = 0.0d0 
vtx_surfseg_new(:) = 0
do vv=1,volume_mesh%nvtx
    if (vtx_map(vv) .GT. 0) then 
        vertices_new(vtx_map(vv),:) = volume_mesh%vertices(vv,:)
        vtx_surfseg_new(vtx_map(vv)) = volume_mesh%vtx_surfseg(vv)
    end if 
end do 

!Map updated edges 
do ee=1,Nenew
    edge_new(ee,1) = vtx_map(edge_new(ee,1)) 
    edge_new(ee,2) = vtx_map(edge_new(ee,2)) 
end do 

!Store new items in the volume mesh structure 
deallocate(volume_mesh%edge)
deallocate(volume_mesh%vertices)
deallocate(volume_mesh%vtx_surfseg)
allocate(volume_mesh%edge(Nenew,4))
volume_mesh%edge(:,:) = edge_new(:,:)
allocate(volume_mesh%vertices(Nvnew,2))
volume_mesh%vertices(:,:) = vertices_new(:,:)
allocate(volume_mesh%vtx_surfseg(Nvnew))
volume_mesh%vtx_surfseg(:) = vtx_surfseg_new(:)
volume_mesh%nedge = Nenew
volume_mesh%nvtx = Nvnew

!Display
if (cm2dopt%dispt == 1) then
    write(*,'(A,I0,A)') '    {collapsed ',nintv2,' internal valence two vertices}'
end if 
return 
end subroutine clean_internal_vlnc2_vertices




!Subroutine divide cells to ensure only tri and quad cells exist in the mesh ===========================
subroutine split_vlnc_gt4cells(volume_mesh)
implicit none 

!Variables - Import
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ee,cc
integer(in) :: c3,c4,maxcedge,ncsplit,Nvnew,Ncnew,Nenew,Vins,Cins,Eins,etgt,ctgt,v1,v2,en,ep
integer(in) :: nedge_cell(volume_mesh%ncell),cell2split(volume_mesh%ncell),cell_vtx_idx(volume_mesh%ncell)
integer(in) :: V2E(volume_mesh%nvtx,2)
integer(in), dimension(:), allocatable :: cedge_loop,cvtx_loop,cell_new,clevtemp
integer(in), dimension(:,:), allocatable :: cell2edge,edge_new
real(dp) :: EWsum,Ledge
real(dp) :: vmtemp(2)
real(dp), dimension(:,:), allocatable :: vertices_new,edge_midpoint_new

!Build list of edges in each cell 
nedge_cell(:) = 0 
do ee=1,volume_mesh%nedge
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    if (c3 .GT. 0) then 
        nedge_cell(c3) = nedge_cell(c3) + 1
    end if 
    if (c4 .GT. 0) then 
        nedge_cell(c4) = nedge_cell(c4) + 1
    end if 
end do 
maxcedge = maxval(nedge_cell)
allocate(cell2edge(volume_mesh%ncell,maxcedge))
cell2edge(:,:) = 0 
nedge_cell(:) = 0 
do ee=1,volume_mesh%nedge
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    if (c3 .GT. 0) then 
        nedge_cell(c3) = nedge_cell(c3) + 1
        cell2edge(c3,nedge_cell(c3)) = ee 
    end if 
    if (c4 .GT. 0) then 
        nedge_cell(c4) = nedge_cell(c4) + 1
        cell2edge(c4,nedge_cell(c4)) = ee 
    end if 
end do 

!Identify any cells that must be split 
ncsplit = 0 
cell2split(:) = 0 
do cc=1,volume_mesh%ncell
    if (nedge_cell(cc) .GT. 4) then 
        ncsplit = ncsplit + 1
        cell2split(cc) = 1
    end if 
end do 

!Set new vertex count 
Nvnew = volume_mesh%nvtx + ncsplit
allocate(vertices_new(Nvnew,2))
vertices_new(:,:) = 0.0d0 
vertices_new(1:volume_mesh%nvtx,:) = volume_mesh%vertices(:,:)

!Set new cell count 
Ncnew = volume_mesh%ncell 
do cc=1,volume_mesh%ncell
    if (cell2split(cc) == 1) then 
        Ncnew = Ncnew + nedge_cell(cc)
    end if 
end do 
allocate(clevtemp(Ncnew))
clevtemp(:) = 0
clevtemp(1:volume_mesh%ncell) = volume_mesh%cell_level(1:volume_mesh%ncell)

!Set new edge count 
Nenew = volume_mesh%nedge
do cc=1,volume_mesh%ncell
    if (cell2split(cc) == 1) then 
        Nenew = Nenew + nedge_cell(cc) 
    end if
end do 
allocate(edge_new(Nenew,4))
allocate(edge_midpoint_new(Nenew,2))
edge_new(:,:) = 0 
edge_new(1:volume_mesh%nedge,:) = volume_mesh%edge(:,:)
edge_midpoint_new(:,:) = 0.0d0 
edge_midpoint_new(1:volume_mesh%nedge,:) = volume_mesh%edge_midpoint(:,:)

!Index split vertices 
Vins = volume_mesh%nvtx 
cell_vtx_idx(:) = 0 
do cc=1,volume_mesh%ncell
    if (cell2split(cc) == 1) then 
        Vins = Vins + 1
        cell_vtx_idx(cc) = Vins
    end if
end do 

!Construct new vertices 
do cc=1,volume_mesh%ncell
    if (cell2split(cc) == 1) then 
        EWsum = 0.0d0 
        vmtemp(:) = 0.0d0 
        do ee=1,nedge_cell(cc)
            etgt = cell2edge(cc,ee)
            v1 = volume_mesh%edge(etgt,1)
            v2 = volume_mesh%edge(etgt,2)
            Ledge = norm2(volume_mesh%vertices(v2,:) - volume_mesh%vertices(v1,:))
            EWsum = EWsum + Ledge 
            vmtemp(1) = vmtemp(1) + 0.5d0*Ledge*(volume_mesh%vertices(v2,1) + volume_mesh%vertices(v1,1))
            vmtemp(2) = vmtemp(2) + 0.5d0*Ledge*(volume_mesh%vertices(v2,2) + volume_mesh%vertices(v1,2))
        end do 
        vmtemp(:) = vmtemp(:)/EWsum
        vertices_new(cell_vtx_idx(cc),:) = vmtemp(:)
    end if 
end do 

!Subdivide targeted cells 
Eins = volume_mesh%nedge
Cins = volume_mesh%ncell
V2E(:,:) = 0 
allocate(cedge_loop(maxcedge)) 
allocate(cvtx_loop(maxcedge))
allocate(cell_new(maxcedge))
do cc=1,volume_mesh%ncell
    if (cell2split(cc) == 1) then 

        !Build edges and vertices for this cell 
        call build_cell_edges_vertices(cedge_loop,cvtx_loop,volume_mesh,nedge_cell,cell2edge,V2E,cc)

        !Add new cells on each looped edge 
        cell_new(:) = 0 
        do ee=1,nedge_cell(cc)
            Cins = Cins + 1
            cell_new(ee) = Cins 
            clevtemp(Cins) = volume_mesh%cell_level(cc)
        end do 

        !Build new edges creating the new cells 
        do ee=1,nedge_cell(cc)
            
            !Increment edge 
            Eins = Eins + 1
            
            !Left and right cells 
            en = ee 
            ep = modulo(ee-2,nedge_cell(cc)) + 1
            c3 = cell_new(en)
            c4 = cell_new(ep)

            !Add edge 
            edge_new(Eins,1) = cell_vtx_idx(cc)
            edge_new(Eins,2) = cvtx_loop(en)
            edge_new(Eins,3) = c3 
            edge_new(Eins,4) = c4 

            !Update cell state on edge en to it's new cell 
            etgt = cedge_loop(ep)
            ctgt = cell_new(ep)
            if (edge_new(etgt,3) == cc) then 
                edge_new(etgt,3) = ctgt
            end if 
            if (edge_new(etgt,4) == cc) then 
                edge_new(etgt,4) = ctgt
            end if 
        end do 

        !Reset V2E
        do ee=1,nedge_cell(cc)
            etgt = cell2edge(cc,ee)
            v1 = volume_mesh%edge(etgt,1)
            v2 = volume_mesh%edge(etgt,2)
            V2E(v1,:) = 0 
            V2E(v2,:) = 0
        end do 
    end if 
end do 

!Store new items in the volume mesh structure 
deallocate(volume_mesh%edge)
deallocate(volume_mesh%vertices)
deallocate(volume_mesh%cell_level)
deallocate(volume_mesh%edge_midpoint)
allocate(volume_mesh%edge(Nenew,4))
volume_mesh%edge(:,:) = edge_new(:,:)
allocate(volume_mesh%edge_midpoint(Nenew,2))
volume_mesh%edge_midpoint(:,:) = edge_midpoint_new(:,:)
allocate(volume_mesh%vertices(Nvnew,2))
volume_mesh%vertices(:,:) = vertices_new(:,:)
allocate(volume_mesh%cell_level(Ncnew))
volume_mesh%cell_level(:) = clevtemp(:)
volume_mesh%nedge = Nenew
volume_mesh%nvtx = Nvnew
volume_mesh%ncell = Ncnew
return 
end subroutine split_vlnc_gt4cells




!Subroutine to structure mesh for SU2_dual output format ===========================
subroutine construct_dual_mesh(volume_mesh)
implicit none 

!Variables - Import
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: vv,ee,cc 
integer(in) :: v1,v2,c3,c4,e1,e2,bc1,bc2,bcnew,Ncnew,Nvnew,Nenew,maxcedge
integer(in) :: cell_vtxidx(volume_mesh%ncell),edge_vtxidx(volume_mesh%nedge),vtx_cidx(volume_mesh%nvtx),vtx_vidx(volume_mesh%nvtx)
integer(in) :: vtx_bndry(volume_mesh%nvtx),vtx_onesurf(volume_mesh%nvtx),nedge_cell(volume_mesh%ncell)
integer(in) :: edge_new(2*volume_mesh%nedge,4),V2E(volume_mesh%nvtx,2)
integer(in), dimension(:,:), allocatable :: cell2edge
real(dp) :: emx,emy,ledge
real(dp) :: cvtx_new_weight(volume_mesh%ncell)
real(dp) :: cvtx_new(volume_mesh%ncell,2),evtx_new(volume_mesh%nedge,2)
real(dp), dimension(:,:), allocatable :: vertices_new

!Build list of edges in each cell 
nedge_cell(:) = 0 
do ee=1,volume_mesh%nedge
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    if (c3 .GT. 0) then 
        nedge_cell(c3) = nedge_cell(c3) + 1
    end if 
    if (c4 .GT. 0) then 
        nedge_cell(c4) = nedge_cell(c4) + 1
    end if 
end do 
maxcedge = maxval(nedge_cell)
allocate(cell2edge(volume_mesh%ncell,maxcedge))
cell2edge(:,:) = 0 
nedge_cell(:) = 0 
do ee=1,volume_mesh%nedge
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    if (c3 .GT. 0) then 
        nedge_cell(c3) = nedge_cell(c3) + 1
        cell2edge(c3,nedge_cell(c3)) = ee 
    end if 
    if (c4 .GT. 0) then 
        nedge_cell(c4) = nedge_cell(c4) + 1
        cell2edge(c4,nedge_cell(c4)) = ee 
    end if 
end do 

!Tag surface vertices 
vtx_bndry(:) = 0
do ee=1,volume_mesh%nedge !all surface and boundary edges are ordered 
    v1 = volume_mesh%edge(ee,1)
    v2 = volume_mesh%edge(ee,2)
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    if ((c3 .LT. 0) .OR. (c4 .LT. 0)) then
        vtx_bndry(v1) = 1
        vtx_bndry(v2) = 1
    end if
end do 

!Tag vertices on edges attached to surface geometry 
vtx_onesurf(:) = 0
do ee=1,volume_mesh%nedge !all surface and boundary edges are ordered 
    v1 = volume_mesh%edge(ee,1)
    v2 = volume_mesh%edge(ee,2)
    if ((vtx_bndry(v1) == 1) .OR. (vtx_bndry(v2) == 1)) then
        vtx_onesurf(v1) = 1
        vtx_onesurf(v2) = 1
    end if 
end do 

!Build V2E for boundary condition vertices and edges
V2E(:,:) = 0.0d0 
do ee=1,volume_mesh%nedge !all surface and boundary edges are ordered 
    v1 = volume_mesh%edge(ee,1)
    v2 = volume_mesh%edge(ee,2)
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    if ((c3 .LT. 0) .OR. (c4 .LT. 0)) then
        V2E(v1,2) = ee
        V2E(v2,1) = ee 
    end if 
end do 

!Index new cells at mesh vertex locations 
Ncnew = 0 
vtx_cidx(:) = 0 
do vv=1,volume_mesh%nvtx
    Ncnew = Ncnew + 1
    vtx_cidx(vv) = Ncnew
end do 

!Index new vertices at vertex, cell and edge locations 
Nvnew = 0 
cell_vtxidx(:) = 0 
edge_vtxidx(:) = 0 
vtx_vidx(:) = 0 
do cc=1,volume_mesh%ncell
    Nvnew = Nvnew + 1
    cell_vtxidx(cc) = Nvnew
end do 
do ee=1,volume_mesh%nedge
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    if ((c3 .LT. 0) .OR. (c4 .LT. 0)) then
        Nvnew = Nvnew + 1
        edge_vtxidx(ee) = Nvnew
    end if
end do 
do vv=1,volume_mesh%nvtx 
    if (vtx_bndry(vv) == 1) then 

        !Adjacent edges 
        e1 = V2E(vv,1)
        e2 = V2E(vv,2)

        !Adjacent boundary conditions 
        bc1 = 0 
        if (volume_mesh%edge(e1,3) .LT. 0) then 
            bc1 = volume_mesh%edge(e1,3)
        elseif (volume_mesh%edge(e1,4) .LT. 0) then 
            bc1 = volume_mesh%edge(e1,4)
        end if 
        bc2 = 0 
        if (volume_mesh%edge(e2,3) .LT. 0) then 
            bc2 = volume_mesh%edge(e2,3)
        elseif (volume_mesh%edge(e2,4) .LT. 0) then 
            bc2 = volume_mesh%edge(e2,4)
        end if 

        !Retain if differing boundary conditions 
        if (bc1 .NE. bc2) then 
            Nvnew = Nvnew + 1
            vtx_vidx(vv) = Nvnew
        end if 
    end if
end do 

!Build new vertices in current cells and edge midpoints 
cvtx_new(:,:) = 0.0d0 
evtx_new(:,:) = 0.0d0 
cvtx_new_weight(:) = 0.0d0 
do ee=1,volume_mesh%nedge
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    emx = 0.5d0*(volume_mesh%vertices(volume_mesh%edge(ee,1),1) + volume_mesh%vertices(volume_mesh%edge(ee,2),1))
    emy = 0.5d0*(volume_mesh%vertices(volume_mesh%edge(ee,1),2) + volume_mesh%vertices(volume_mesh%edge(ee,2),2))
    ledge = norm2(volume_mesh%vertices(volume_mesh%edge(ee,2),:) - volume_mesh%vertices(volume_mesh%edge(ee,1),:))
    evtx_new(ee,1) = emx 
    evtx_new(ee,2) = emy 
    if (c3 .GT. 0) then 
        cvtx_new(c3,1) = cvtx_new(c3,1) + emx*ledge
        cvtx_new(c3,2) = cvtx_new(c3,2) + emy*ledge
        cvtx_new_weight(c3) = cvtx_new_weight(c3) + ledge
    end if 
    if (c4 .GT. 0) then 
        cvtx_new(c4,1) = cvtx_new(c4,1) + emx*ledge
        cvtx_new(c4,2) = cvtx_new(c4,2) + emy*ledge
        cvtx_new_weight(c4) = cvtx_new_weight(c4) + ledge
    end if 
end do 
cvtx_new(:,1) = cvtx_new(:,1)/cvtx_new_weight(:)
cvtx_new(:,2) = cvtx_new(:,2)/cvtx_new_weight(:)

!Build new vertex list 
allocate(vertices_new(Nvnew,2))
vertices_new(:,:) = 0.0d0 
do cc=1,volume_mesh%ncell
    vertices_new(cell_vtxidx(cc),:) = cvtx_new(cc,:)
end do 
do ee=1,volume_mesh%nedge
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    if ((c3 .LT. 0) .OR. (c4 .LT. 0)) then
        vertices_new(edge_vtxidx(ee),:) = evtx_new(ee,:)
    end if
end do 
do vv=1,volume_mesh%nvtx 
    if (vtx_vidx(vv) .NE. 0) then 
        vertices_new(vtx_vidx(vv),:) = volume_mesh%vertices(vv,:)
    end if
end do 

!Build new edges connecting cell and edge midpoint vertices
Nenew = 0
edge_new(:,:) = 0 
do ee=1,volume_mesh%nedge
    v1 = volume_mesh%edge(ee,1)
    v2 = volume_mesh%edge(ee,2)
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    if ((c3 .GT. 0) .AND. (c4 .GT. 0)) then
        Nenew = Nenew + 1
        edge_new(Nenew,1) = cell_vtxidx(c3)
        edge_new(Nenew,2) = cell_vtxidx(c4)
        edge_new(Nenew,3) = vtx_cidx(v2)
        edge_new(Nenew,4) = vtx_cidx(v1)
    elseif (c3 .LT. 0) then
        Nenew = Nenew + 1
        edge_new(Nenew,1) = cell_vtxidx(c4)
        edge_new(Nenew,2) = edge_vtxidx(ee)
        edge_new(Nenew,3) = vtx_cidx(v1)
        edge_new(Nenew,4) = vtx_cidx(v2)
    elseif (c4 .LT. 0) then
        Nenew = Nenew + 1
        edge_new(Nenew,1) = cell_vtxidx(c3)
        edge_new(Nenew,2) = edge_vtxidx(ee)
        edge_new(Nenew,3) = vtx_cidx(v1)
        edge_new(Nenew,4) = vtx_cidx(v2)
    end if 
end do 

!Build edges connecting edges either side of boundary condition vertices 
do vv=1,volume_mesh%nvtx 
    if (vtx_bndry(vv) == 1) then 

        !Adjacent edges 
        e1 = V2E(vv,1)
        e2 = V2E(vv,2)

        !Adjacent boundary conditions 
        bc1 = 0 
        if (volume_mesh%edge(e1,3) .LT. 0) then 
            bc1 = volume_mesh%edge(e1,3)
        elseif (volume_mesh%edge(e1,4) .LT. 0) then 
            bc1 = volume_mesh%edge(e1,4)
        end if 
        bc2 = 0 
        if (volume_mesh%edge(e2,3) .LT. 0) then 
            bc2 = volume_mesh%edge(e2,3)
        elseif (volume_mesh%edge(e2,4) .LT. 0) then 
            bc2 = volume_mesh%edge(e2,4)
        end if 

        !Edge construction 
        if (bc1 == bc2) then !build one edge
            bcnew = bc1 
            Nenew = Nenew + 1
            edge_new(Nenew,1) = edge_vtxidx(e1)
            edge_new(Nenew,2) = edge_vtxidx(e2)
            edge_new(Nenew,3) = bcnew
            edge_new(Nenew,4) = vtx_cidx(vv)
        else !build two edges 
            Nenew = Nenew + 1
            edge_new(Nenew,1) = edge_vtxidx(e1)
            edge_new(Nenew,2) = vtx_vidx(vv)
            edge_new(Nenew,3) = bc1
            edge_new(Nenew,4) = vtx_cidx(vv)
            Nenew = Nenew + 1
            edge_new(Nenew,1) = vtx_vidx(vv)
            edge_new(Nenew,2) = edge_vtxidx(e2)
            edge_new(Nenew,3) = bc2
            edge_new(Nenew,4) = vtx_cidx(vv)
        end if 
    end if 
end do 

!Store new items in the volume mesh structure 
deallocate(volume_mesh%edge)
deallocate(volume_mesh%vertices)
deallocate(volume_mesh%cell_level)
allocate(volume_mesh%edge(Nenew,4))
volume_mesh%edge(:,:) = edge_new(1:Nenew,:)
allocate(volume_mesh%vertices(Nvnew,2))
volume_mesh%vertices(:,:) = vertices_new(:,:)
allocate(volume_mesh%cell_level(Ncnew))
volume_mesh%cell_level(:) = 0
volume_mesh%nedge = Nenew
volume_mesh%nvtx = Nvnew
volume_mesh%ncell = Ncnew
! open(11,file='io/vtxtest.dat')
! do vv=1,Nvnew
!     write(11,*) vertices_new(vv,:)
! end do 
! close(11)
return 
end subroutine construct_dual_mesh




!Deform mesh to surface subroutine ===========================
subroutine deform_mesh2surface(volume_mesh,surface_mesh,surface_adtree)
implicit none 

!Variables - Import
type(vol_mesh_data) :: volume_mesh
type(surface_data) :: surface_mesh
type(tree_data) :: surface_adtree

!Variables - Local 
integer(in) :: ee,vv,aa,nn,kk,ff
integer(in) :: Nsurf_edge,etgt,maxvlnc,v1,v2,nselected
integer(in) :: vtx_surface(volume_mesh%nvtx)
integer(in) :: valence(volume_mesh%nvtx),node_select(surface_adtree%nnode)
integer(in), dimension(:), allocatable :: vmsurf_edges
integer(in), dimension(:,:), allocatable :: v2v
real(dp) :: dx,dy,nnorm,xmin,ymin,xmax,ymax,rpad,rsup,normN,normC,adj_weight,adj_weightT,dref
real(dp) :: zxmin,zxmax,zymin,zymax,zzmin,zzmax
real(dp) :: vl1(2),vl2(2),np1(2),np2(2),vi(2),vic(2)
real(dp) :: vtx_normal(volume_mesh%nvtx,2),vtx_deformation(volume_mesh%nvtx,2),vtx_deformationN(volume_mesh%nvtx,2)

!Build mesh valence
valence(:) = 0 
do ee=1,volume_mesh%nedge
    valence(volume_mesh%edge(ee,1:2)) = valence(volume_mesh%edge(ee,1:2)) + 1
end do 
maxvlnc = maxval(valence)

!Build mesh v2v
allocate(v2v(volume_mesh%nvtx,maxvlnc))
v2v(:,:) = 0
do ee=1,volume_mesh%nedge
    v1 = volume_mesh%edge(ee,1)
    v2 = volume_mesh%edge(ee,2)
    do aa=1,maxvlnc
        if (v2v(v1,aa) == 0) then 
            v2v(v1,aa) = v2 
            exit
        end if 
    end do 
    do aa=1,maxvlnc
        if (v2v(v2,aa) == 0) then 
            v2v(v2,aa) = v1
            exit
        end if 
    end do 
end do 

!Identify surface edges in the volume mesh
Nsurf_edge = 0 
do ee=1,volume_mesh%nedge
    if ((volume_mesh%edge(ee,3) == -1) .OR. (volume_mesh%edge(ee,4) == -1)) then 
        Nsurf_edge = Nsurf_edge + 1
    end if 
end do 
allocate(vmsurf_edges(Nsurf_edge))
Nsurf_edge = 0 
do ee=1,volume_mesh%nedge
    if ((volume_mesh%edge(ee,3) == -1) .OR. (volume_mesh%edge(ee,4) == -1)) then 
        Nsurf_edge = Nsurf_edge + 1
        vmsurf_edges(Nsurf_edge) = ee 
    end if 
end do 

!Identify surface vertices
vtx_surface(:) = 0 
do ee=1,Nsurf_edge
    etgt = vmsurf_edges(ee)
    vtx_surface(volume_mesh%edge(etgt,1:2)) = 1 
end do 

!Find surface normal directions in the volume mesh surface at each vertex
vtx_normal(:,:) = 0.0d0 
do ee=1,Nsurf_edge
    etgt = vmsurf_edges(ee)
    dx = volume_mesh%vertices(volume_mesh%edge(etgt,2),1) - volume_mesh%vertices(volume_mesh%edge(etgt,1),1)
    dy = volume_mesh%vertices(volume_mesh%edge(etgt,2),2) - volume_mesh%vertices(volume_mesh%edge(etgt,1),2)
    vtx_normal(volume_mesh%edge(etgt,1),1) = vtx_normal(volume_mesh%edge(etgt,1),1) + dy
    vtx_normal(volume_mesh%edge(etgt,1),2) = vtx_normal(volume_mesh%edge(etgt,1),2) - dx
    vtx_normal(volume_mesh%edge(etgt,2),1) = vtx_normal(volume_mesh%edge(etgt,2),1) + dy
    vtx_normal(volume_mesh%edge(etgt,2),2) = vtx_normal(volume_mesh%edge(etgt,2),2) - dx
end do 
do vv=1,volume_mesh%nvtx
    if (vtx_surface(vv) == 1) then 
        nnorm = norm2(vtx_normal(vv,:))
        if (nnorm .NE. 0.0d0) then 
            vtx_normal(vv,:) = vtx_normal(vv,:)/nnorm
        end if 
    end if 
end do 

!Find surface mesh intersection and hence deformation for each surface vertex
zzmin = 0.0d0 
zzmax = 0.0d0 
nselected = 0 
node_select(:) = 0 
vtx_deformation(:,:) = 0.0d0 
do vv=1,volume_mesh%nvtx
    if (vtx_surface(vv) == 1) then 

        !Build search radius 
        xmax = max(maxval(volume_mesh%vertices(v2v(vv,1:valence(vv)),1)),volume_mesh%vertices(vv,1))
        ymax = max(maxval(volume_mesh%vertices(v2v(vv,1:valence(vv)),2)),volume_mesh%vertices(vv,2))
        xmin = min(minval(volume_mesh%vertices(v2v(vv,1:valence(vv)),1)),volume_mesh%vertices(vv,1))
        ymin = min(minval(volume_mesh%vertices(v2v(vv,1:valence(vv)),2)),volume_mesh%vertices(vv,2))
        rpad = sqrt((xmax - xmin)**2 + (ymax - ymin)**2)

        !Build padded zone 
        zxmin = volume_mesh%vertices(vv,1) - rpad
        zxmax = volume_mesh%vertices(vv,1) + rpad
        zymin = volume_mesh%vertices(vv,2) - rpad
        zymax = volume_mesh%vertices(vv,2) + rpad

        !Build normal search line 
        np1(:) = volume_mesh%vertices(vv,:) - rpad*vtx_normal(vv,:)
        np2(:) = volume_mesh%vertices(vv,:) + rpad*vtx_normal(vv,:)

        !Search for surface mesh items 
        call search_ADtree(nselected,node_select,surface_adtree,zxmin,zxmax,zymin,zymax,zzmin,zzmax)

        !If any found then search for closest intersection 
        if (nselected .NE. 0) then 
            vic(:) = ieee_value(1.0d0,IEEE_POSITIVE_INF)
            do nn=1,nselected
                do kk=1,surface_adtree%tree(node_select(nn))%nentry
            
                    !Verticies of the ends of the surface segment 
                    vl1(:) = surface_mesh%vertices(surface_mesh%faces(surface_adtree%tree(node_select(nn))%entry(kk),1),:)
                    vl2(:) = surface_mesh%vertices(surface_mesh%faces(surface_adtree%tree(node_select(nn))%entry(kk),2),:)

                    !Check for intersection 
                    vi = line_line_intersection_loc_inl1(vl1,vl2,np1,np2)

                    !If valid 
                    if (.NOT.isnan(vi(1))) then 
                        normN = sqrt((vi(1) - volume_mesh%vertices(vv,1))**2 + (vi(2) - volume_mesh%vertices(vv,2))**2)
                        normC = sqrt((vic(1) - volume_mesh%vertices(vv,1))**2 + (vic(2) - volume_mesh%vertices(vv,2))**2)
                        if (normN .LT. normC) then 
                            vic = vi
                        end if 
                    end if 
                end do 
            end do 
            if (.NOT.isnan(vic(1))) then 
                vtx_deformation(vv,:) = vic(:) - volume_mesh%vertices(vv,:)
            end if 
        end if 
    end if 
end do 

!Flood deformations through mesh 
do ff=1,20
    do vv=1,volume_mesh%nvtx
        if (vtx_surface(vv) == 1) then !fixed surface deformation
            vtx_deformationN(vv,:) = vtx_deformation(vv,:)
        else !interpolated volume deformation 

            !Build support radius 
            xmax = max(maxval(volume_mesh%vertices(v2v(vv,1:valence(vv)),1)),volume_mesh%vertices(vv,1))
            ymax = max(maxval(volume_mesh%vertices(v2v(vv,1:valence(vv)),2)),volume_mesh%vertices(vv,2))
            xmin = min(minval(volume_mesh%vertices(v2v(vv,1:valence(vv)),1)),volume_mesh%vertices(vv,1))
            ymin = min(minval(volume_mesh%vertices(v2v(vv,1:valence(vv)),2)),volume_mesh%vertices(vv,2))
            rsup = 1.1d0*sqrt((xmax - xmin)**2 + (ymax - ymin)**2)

            !Accumulate deformation 
            adj_weightT = 0.0d0 
            vtx_deformationN(vv,:) = 0.0d0 
            do aa=1,valence(vv)
                dref = norm2(volume_mesh%vertices(v2v(vv,aa),:) + vtx_deformation(v2v(vv,aa),:) - &
                volume_mesh%vertices(vv,:) - vtx_deformation(vv,:))
                adj_weight = wendlandc2(dref,rsup)
                vtx_deformationN(vv,:) = vtx_deformationN(vv,:) + adj_weight*vtx_deformation(v2v(vv,aa),:)
                adj_weightT = adj_weightT + adj_weight
            end do 
            if (adj_weightT .NE. 0.0d0) then 
                vtx_deformationN(vv,:) = vtx_deformationN(vv,:)/adj_weightT
            end if 
        end if 
    end do 
    vtx_deformation(:,:) = vtx_deformationN(:,:)
end do 

!Snap mesh to the geometry surface 
volume_mesh%vertices(:,:) = volume_mesh%vertices(:,:) + vtx_deformation(:,:)
return 
end subroutine deform_mesh2surface




!Build mesh cells ===========================
subroutine build_mesh_cells(volume_mesh)
implicit none 

!Variables - Import
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ee,cc
integer(in) :: c3,c4,maxcedge,Cins,Eins,etgt,v1,v2
integer(in) :: nedge_cell(volume_mesh%ncell)
integer(in) :: V2E(volume_mesh%nvtx,2)
integer(in), dimension(:), allocatable :: cedge_loop,cvtx_loop
integer(in), dimension(:,:), allocatable :: cell2edge

!Build list of edges in each cell 
nedge_cell(:) = 0 
do ee=1,volume_mesh%nedge
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    if (c3 .GT. 0) then 
        nedge_cell(c3) = nedge_cell(c3) + 1
    end if 
    if (c4 .GT. 0) then 
        nedge_cell(c4) = nedge_cell(c4) + 1
    end if 
end do 
maxcedge = maxval(nedge_cell)
allocate(cell2edge(volume_mesh%ncell,maxcedge))
cell2edge(:,:) = 0 
nedge_cell(:) = 0 
do ee=1,volume_mesh%nedge
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    if (c3 .GT. 0) then 
        nedge_cell(c3) = nedge_cell(c3) + 1
        cell2edge(c3,nedge_cell(c3)) = ee 
    end if 
    if (c4 .GT. 0) then 
        nedge_cell(c4) = nedge_cell(c4) + 1
        cell2edge(c4,nedge_cell(c4)) = ee 
    end if 
end do 

!Allocate cell structure 
allocate(volume_mesh%cells(volume_mesh%ncell))

!Construct cells 
Eins = volume_mesh%nedge
Cins = volume_mesh%ncell
V2E(:,:) = 0 
allocate(cedge_loop(maxcedge)) 
allocate(cvtx_loop(maxcedge))
do cc=1,volume_mesh%ncell

    !Build edges and vertices for this cell 
    call build_cell_edges_vertices(cedge_loop,cvtx_loop,volume_mesh,nedge_cell,cell2edge,V2E,cc)

    !Store data
    volume_mesh%cells(cc)%nvtx = nedge_cell(cc)
    volume_mesh%cells(cc)%nedge = nedge_cell(cc)
    allocate(volume_mesh%cells(cc)%vertices(nedge_cell(cc)))
    allocate(volume_mesh%cells(cc)%edges(nedge_cell(cc)))
    volume_mesh%cells(cc)%vertices(:) = cvtx_loop(1:nedge_cell(cc))
    volume_mesh%cells(cc)%edges(:) = cedge_loop(1:nedge_cell(cc))

    !Reset V2E
    do ee=1,nedge_cell(cc)
        etgt = cell2edge(cc,ee)
        v1 = volume_mesh%edge(etgt,1)
        v2 = volume_mesh%edge(etgt,2)
        V2E(v1,:) = 0 
        V2E(v2,:) = 0
    end do 
end do 
return 
end subroutine build_mesh_cells




!Subroutine to build mesh cell edge and vertex loops ===========================
subroutine build_cell_edges_vertices(cedge_loop,cvtx_loop,volume_mesh,nedge_cell,cell2edge,V2E,ctgt)
implicit none 

!Variables - Import
integer(in) :: ctgt 
integer(in), dimension(:) :: cedge_loop,cvtx_loop,nedge_cell
integer(in), dimension(:,:) :: cell2edge,V2E
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ee,aa,vv,ii
integer(in) :: etgt,eadj,eadjv,vtgt,v1,v2,en,ep
integer(in) :: cedge_tag(size(cedge_loop,1)),loop_flip(size(cedge_loop,1))
real(dp) :: Acell

!Build V2E for existing edges on this cell 
do ee=1,nedge_cell(ctgt)
    etgt = cell2edge(ctgt,ee)
    v1 = volume_mesh%edge(etgt,1)
    v2 = volume_mesh%edge(etgt,2)
    do aa=1,2
        if (V2E(v1,aa) == 0) then 
            V2E(v1,aa) = etgt 
            exit 
        end if 
    end do 
    do aa=1,2
        if (V2E(v2,aa) == 0) then 
            V2E(v2,aa) = etgt 
            exit 
        end if 
    end do 
end do 

!Initialise edge loop 
cedge_tag(:) = 0 
cedge_loop(:) = 0 
cedge_tag(1) = 1
cedge_loop(1) = cell2edge(ctgt,1)
etgt = cell2edge(ctgt,1)

!Order edges around this cell 
do ee=2,nedge_cell(ctgt)

    !Find untagged edge adjacent to the current edge 
    eadjv = 0 
    do vv=1,2
        vtgt = volume_mesh%edge(etgt,vv)
        do aa=1,2
            eadj = V2E(vtgt,aa)
            if (eadj .GT. 0) then 
                do ii=1,nedge_cell(ctgt)
                    if ((cell2edge(ctgt,ii) == eadj) .AND. (cedge_tag(ii) == 0)) then 
                        eadjv = eadj 
                        cedge_tag(ii) = 1
                        exit 
                    end if 
                end do 
                if (eadjv .NE. 0) then 
                    exit 
                end if 
            end if 
        end do 
        if (eadjv .NE. 0) then 
            exit 
        end if 
    end do 

    !Exit if no edge found as loop is complete 
    if (eadjv == 0) then 
        exit 
    end if 

    !Add new edge
    cedge_loop(ee) = eadjv

    !Update current edge 
    etgt = eadjv 
end do 

!Build vertex loop 
cvtx_loop(:) = 0
do ee=1,nedge_cell(ctgt)

    !Next and previous edge
    en = ee 
    ep = modulo(ee-2,nedge_cell(ctgt)) + 1
    en = cedge_loop(en)
    ep = cedge_loop(ep)
    
    !Find shared vertex
    vtgt = 0 
    do vv=1,2
        v1 = volume_mesh%edge(en,vv)
        if ((volume_mesh%edge(ep,1) == v1) .OR. (volume_mesh%edge(ep,2) == v1)) then 
            vtgt = v1 
            exit 
        end if 
    end do 

    !Store
    cvtx_loop(ee) = vtgt
end do  

!Check orientation 
Acell = 0.0d0 
do ee=1,nedge_cell(ctgt)
    v1 = ee 
    v2 = modulo(ee,nedge_cell(ctgt)) + 1
    v1 = cvtx_loop(v1)
    v2 = cvtx_loop(v2)
    Acell = Acell + Asegment(volume_mesh%vertices(v1,:),volume_mesh%vertices(v2,:))
end do 

!Flip if required 
if (Acell .LT. 0.0d0) then 
    loop_flip(:) = cedge_loop(:)
    do ee=1,nedge_cell(ctgt)
        cedge_loop(ee) = loop_flip(nedge_cell(ctgt)-ee+1)
    end do 
    cvtx_loop(:) = 0
    do ee=1,nedge_cell(ctgt)

        !Next and previous edge
        en = ee 
        ep = modulo(ee-2,nedge_cell(ctgt)) + 1
        en = cedge_loop(en)
        ep = cedge_loop(ep)
        
        !Find shared vertex
        vtgt = 0 
        do vv=1,2
            v1 = volume_mesh%edge(en,vv)
            if ((volume_mesh%edge(ep,1) == v1) .OR. (volume_mesh%edge(ep,2) == v1)) then 
                vtgt = v1 
                exit 
            end if 
        end do 

        !Store
        cvtx_loop(ee) = vtgt
    end do  
end if 
return 
end subroutine build_cell_edges_vertices


end module cellmesh2d_postprocess_mod