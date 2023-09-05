!cell_mesh2d quadtree module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 0.2.1
!Updated 15-08-2023

!Module
module cellmesh2d_quadtree_mod
use cellmesh2d_adtree_mod
use cellmesh2d_geometry_mod
contains 


!Quadtree mesh refinement subroutine ===========================
subroutine quadtree_mesh_refine(qt_mesh,surface_mesh,surface_adtree,cm2dopt,obj_cx,obj_cy)
implicit none 

!Variables - Import
real(dp) :: obj_cx,obj_cy
type(cm2d_options) :: cm2dopt
type(quadtree_data) :: qt_mesh
type(surface_data) :: surface_mesh
type(tree_data) :: surface_adtree

!Variables - Local 
integer(in) :: rr,cc,ff,ccf,aa 
integer(in) :: cins,vins,nrefineN,nrefine,ncexit,nselected,Nflood,cidx,cadj
integer(in) :: Nflood_refiter(cm2dopt%Nrefine),node_select(surface_adtree%nnode)
integer(in) :: Qrefine(cm2dopt%Ncell_max),QfloodR(cm2dopt%Ncell_max),QfloodRN(cm2dopt%Ncell_max)
integer(in) :: eopposite(4),edges(4,2),ccadj(4,4)
real(dp) :: Nflood_crr,Nflood_dcr,cpadSZ
real(dp) :: zxmin,zxmax,zymin,zymax,zzmin,zzmax

!Define opposite edges
eopposite(1) = 3
eopposite(2) = 4
eopposite(3) = 1
eopposite(4) = 2

!Define cell edges 
edges(1,1) = 1
edges(1,2) = 2
edges(2,1) = 2
edges(2,2) = 3
edges(3,1) = 3
edges(3,2) = 4
edges(4,1) = 4
edges(4,2) = 1

!Define cell child-child adjacency
ccadj(1,1) = -4
ccadj(1,2) = 2
ccadj(1,3) = 4
ccadj(1,4) = -2
ccadj(2,1) = -3
ccadj(2,2) = -1
ccadj(2,3) = 3
ccadj(2,4) = 1
ccadj(3,1) = 2
ccadj(3,2) = -4 
ccadj(3,3) = -2
ccadj(3,4) = 4
ccadj(4,1) = 1
ccadj(4,2) = 3
ccadj(4,3) = -1
ccadj(4,4) = -3

!Initialise refinement arrays
Qrefine(:) = 0
QfloodR(:) = 0
QfloodRN(:) = 0
nrefine = 0 
nrefineN = 0

!Allocate quadtree mesh arrays
allocate(qt_mesh%cell_level(cm2dopt%Ncell_max))
allocate(qt_mesh%cell_adjacent(cm2dopt%Ncell_max,4)) !-2 = far field || -1 = object boundary || else cell index
allocate(qt_mesh%cell_vcnr(cm2dopt%Ncell_max,4))
allocate(qt_mesh%cell_vemid(cm2dopt%Ncell_max,4))
allocate(qt_mesh%cell_child(cm2dopt%Ncell_max,4))
allocate(qt_mesh%vtx(2*cm2dopt%Ncell_max,2))
qt_mesh%cell_level(:) = 0
qt_mesh%cell_adjacent(:,:) = 0
qt_mesh%cell_vcnr(:,:) = 0
qt_mesh%cell_vemid(:,:) = 0
qt_mesh%cell_child(:,:) = 0
qt_mesh%vtx(:,:) = 0.0d0

!Set initial item counter values
cins = 2 !cell
vins = 5 !vertex 

!Construct first cell to define global mesh domain
if (cm2dopt%set_mbounds == 1) then !set custom bounds
    qt_mesh%vtx(1,1) = cm2dopt%mxmin
    qt_mesh%vtx(1,2) = cm2dopt%mymin 
    qt_mesh%vtx(2,1) = cm2dopt%mxmin
    qt_mesh%vtx(2,2) = cm2dopt%mymax 
    qt_mesh%vtx(3,1) = cm2dopt%mxmax
    qt_mesh%vtx(3,2) = cm2dopt%mymax 
    qt_mesh%vtx(4,1) = cm2dopt%mxmax
    qt_mesh%vtx(4,2) = cm2dopt%mymin
    cm2dopt%far_field_bound = max(abs(cm2dopt%mymax - cm2dopt%mymin),abs(cm2dopt%mxmax - cm2dopt%mxmin))
else !set with far field distance
    qt_mesh%vtx(1,1) = obj_cx - cm2dopt%far_field_bound + cm2dopt%om_offset_x
    qt_mesh%vtx(1,2) = obj_cy - cm2dopt%far_field_bound + cm2dopt%om_offset_y
    qt_mesh%vtx(2,1) = obj_cx - cm2dopt%far_field_bound + cm2dopt%om_offset_x
    qt_mesh%vtx(2,2) = obj_cy + cm2dopt%far_field_bound + cm2dopt%om_offset_y
    qt_mesh%vtx(3,1) = obj_cx + cm2dopt%far_field_bound + cm2dopt%om_offset_x
    qt_mesh%vtx(3,2) = obj_cy + cm2dopt%far_field_bound + cm2dopt%om_offset_y
    qt_mesh%vtx(4,1) = obj_cx + cm2dopt%far_field_bound + cm2dopt%om_offset_x
    qt_mesh%vtx(4,2) = obj_cy - cm2dopt%far_field_bound + cm2dopt%om_offset_y
end if 
qt_mesh%cell_vcnr(1,1) = 1
qt_mesh%cell_vcnr(1,2) = 2
qt_mesh%cell_vcnr(1,3) = 3
qt_mesh%cell_vcnr(1,4) = 4
qt_mesh%cell_adjacent(1,1) = -2
qt_mesh%cell_adjacent(1,2) = -2
qt_mesh%cell_adjacent(1,3) = -2
qt_mesh%cell_adjacent(1,4) = -2
qt_mesh%cell_level(1) = 1
zxmin = 0.0d0 
zxmax = 0.0d0 
zymin = 0.0d0 
zymax = 0.0d0 
zzmin = 0.0d0 
zzmax = 0.0d0 

!Construct refinement adjacency flooding value for each iteration 
Nflood_dcr = (real(cm2dopt%Nrefine_flood_i,dp) - real(cm2dopt%Nrefine_flood_f,dp))/real(cm2dopt%Nrefine-1,dp)
Nflood_refiter(:) = 0
Nflood_refiter(1) = cm2dopt%Nrefine_flood_i
Nflood_crr = real(cm2dopt%Nrefine_flood_i,dp)
do rr=2,cm2dopt%Nrefine
    Nflood_crr = Nflood_crr - Nflood_dcr
    Nflood_refiter(rr) = nint(Nflood_crr,in)
end do
Nflood_refiter(cm2dopt%Nrefine) = cm2dopt%Nrefine_flood_f

!Construct quadtree mesh -------------------------------------------------------
ncexit = 0 
nselected = 0 
node_select(:) = 0 
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '--> refining cells: '
    write(*,'(A)') ' level  |  cells  |  Nflood'
    write(*,'(A)') '----------------------------'
    write(*,'(A,I2,A,I8,A,I4)') '   ',0,'   ',1,'     ',Nflood_refiter(1)
end if
do rr=1,cm2dopt%Nrefine + cm2dopt%NrefineB

    !Set padding range base size -> size of cell at this level 
    cpadSZ = 2.0d0*cm2dopt%far_field_bound/(2.0d0**(rr - 1))

    !Identify cells to refine 
    Qrefine(1:cins-1) = 0
    QfloodR(1:nrefine) = 0 
    nrefine = 0
    do cc=1,cins-1
        if (qt_mesh%cell_level(cc) == rr) then !At current maximum level 
            if (qt_mesh%cell_child(cc,1) == 0) then 

                !Intersection bounding box
                zxmin = qt_mesh%vtx(qt_mesh%cell_vcnr(cc,1),1) - cpadSZ*cm2dopt%ADTpadding !tgt bounding box -> xmin
                zxmax = qt_mesh%vtx(qt_mesh%cell_vcnr(cc,3),1) + cpadSZ*cm2dopt%ADTpadding !tgt bounding box -> xmax
                zymin = qt_mesh%vtx(qt_mesh%cell_vcnr(cc,1),2) - cpadSZ*cm2dopt%ADTpadding !tgt bounding box -> ymin
                zymax = qt_mesh%vtx(qt_mesh%cell_vcnr(cc,3),2) + cpadSZ*cm2dopt%ADTpadding !tgt bounding box -> ymax

                !Identify any edge bounding boxes that may overlap the cell 
                call search_ADtree(nselected,node_select,surface_adtree,zxmin,zxmax,zymin,zymax,zzmin,zzmax)

                !Tag cell to refine 
                if (nselected .NE. 0) then 
                    if (rr .LE. cm2dopt%Nrefine) then !Check for and tag geometry surface overlap (normal refinement)
                        Qrefine(cc) = segment_cell_ovlp_bool(qt_mesh,surface_mesh,surface_adtree,node_select,nselected,cc)
                    else !Check for and tag geometry surface overlap with smaller than cell radius of curvature (boosted refinement)
                        Qrefine(cc) = segment_cell_ovlp_wrcurv_bool(qt_mesh,surface_mesh,surface_adtree,node_select,nselected,cc)
                    end if 
                    if (Qrefine(cc) == 1) then 
                        nrefine = nrefine + 1
                        QfloodR(nrefine) = cc
                    end if
                end if 
            end if 
        end if 
    end do 

    !Exit with no refinement requested
    if (nrefine == 0) then 
        if (cm2dopt%dispt == 1) then
            write(*,'(A)')  '   {refinement completed as no cells tagged}'
        end if
        ncexit = 1 
        exit 
    end if 

    !Flood refinement status 
    if (rr .LE. cm2dopt%Nrefine) then 
        Nflood = Nflood_refiter(rr)
    else
        Nflood = cm2dopt%Nrefine_flood_b 
    end if 
    do ff=1,Nflood
        nrefineN = 0 
        do ccf=1,nrefine
            cidx = QfloodR(ccf)
            do aa=1,4 
                cadj = qt_mesh%cell_adjacent(cidx,aa)
                if (cadj .GT. 0) then 
                    if (qt_mesh%cell_child(cadj,1) == 0) then 
                        if (Qrefine(cadj) == 0) then 
                            nrefineN = nrefineN + 1
                            QfloodRN(nrefineN) = cadj
                            Qrefine(cadj) = 1
                        end if 
                    end if 
                end if 
            end do 
        end do 
        QfloodR(1:nrefineN) = QfloodRN(1:nrefineN)
        QfloodRN(1:nrefineN) = 0
        nrefine = nrefineN
    end do 

    !Refine tagged cells 
    do cc=1,cins-1
        if (Qrefine(cc) == 1) then
            call quadtree_cell_refine(qt_mesh,cm2dopt,eopposite,edges,ccadj,vins,cins,cc)
        end if 
    end do 

    !Display
    if (cm2dopt%dispt == 1) then
        write(*,'(A,I2,A,I8,A,I4)') '   ',rr,'   ',cins-1,'     ',Nflood
    end if

    !Display completion of base refinement 
    if (rr == cm2dopt%Nrefine) then 
        if (cm2dopt%dispt == 1) then
            write(*,'(A)')  '   {base refinement completed}'
        end if
    end if 
end do 
if ((rr == cm2dopt%Nrefine + cm2dopt%NrefineB + 1) .AND. (ncexit == 0)) then 
    if (cm2dopt%dispt == 1) then
        write(*,'(A)')  '   {refinement completed at maximum level}'
    end if
end if 

!Assign refined item quantities
qt_mesh%cins = cins 
qt_mesh%vins = vins 
return 
end subroutine quadtree_mesh_refine




!Quadtree refinement subroutine ===========================
subroutine quadtree_cell_refine(qt_mesh,cm2dopt,eopposite,edges,ccadj,vins,cins,ctgt)
implicit none 

!Variables - Import
integer(in) :: vins,cins,ctgt 
integer(in) :: eopposite(4),edges(4,2),ccadj(4,4)
type(quadtree_data) :: qt_mesh
type(cm2d_options) :: cm2dopt

!Variables - Local 
integer(in) :: aa,ccf
integer(in) :: vtxCC,cadj,child_idx

!Build cell centre vertex 
vtxCC = vins 
qt_mesh%vtx(vins,1) = 0.25d0*sum(qt_mesh%vtx(qt_mesh%cell_vcnr(ctgt,:),1))
qt_mesh%vtx(vins,2) = 0.25d0*sum(qt_mesh%vtx(qt_mesh%cell_vcnr(ctgt,:),2))
vins = vins + 1
if (vins .GE. cm2dopt%Ncell_max) then
    cm2dopt%cm2dfailure = 1
    write(*,'(A)') '** maximum vertex count exceeded'
    return 
end if

!Pull adjacent edge midpoint vertices
do aa=1,4 
    cadj = qt_mesh%cell_adjacent(ctgt,aa)
    if (qt_mesh%cell_vemid(ctgt,aa) == 0) then 
        if (cadj .GT. 0) then 
            if (qt_mesh%cell_level(cadj) == qt_mesh%cell_level(ctgt)) then 
                qt_mesh%cell_vemid(ctgt,aa) = qt_mesh%cell_vemid(cadj,eopposite(aa))
            end if 
        end if
    end if 
end do 

!Build required edge midpoint vertices 
do aa=1,4 
    if (qt_mesh%cell_vemid(ctgt,aa) == 0) then 
        qt_mesh%cell_vemid(ctgt,aa) = vins 
        qt_mesh%vtx(vins,1) = 0.5d0*(qt_mesh%vtx(qt_mesh%cell_vcnr(ctgt,edges(aa,1)),1) + &
                                        qt_mesh%vtx(qt_mesh%cell_vcnr(ctgt,edges(aa,2)),1))
        qt_mesh%vtx(vins,2) = 0.5d0*(qt_mesh%vtx(qt_mesh%cell_vcnr(ctgt,edges(aa,1)),2) + &
                                        qt_mesh%vtx(qt_mesh%cell_vcnr(ctgt,edges(aa,2)),2))
        vins = vins + 1
        if (vins .GE. cm2dopt%Ncell_max) then
            cm2dopt%cm2dfailure = 1
            write(*,'(A)') '** maximum vertex count exceeded'
            return 
        end if
    end if 
end do 

!Push edge midpoint vertices 
do aa=1,4 
    cadj = qt_mesh%cell_adjacent(ctgt,aa)
    if (cadj .GT. 0) then 
        if (qt_mesh%cell_level(cadj) == qt_mesh%cell_level(ctgt)) then 
            if (qt_mesh%cell_vemid(cadj,eopposite(aa)) == 0) then 
                qt_mesh%cell_vemid(cadj,eopposite(aa)) = qt_mesh%cell_vemid(ctgt,aa)
            end if
        end if 
    end if 
end do 

!Construct new child cells ----
!child 1
qt_mesh%cell_vcnr(cins,1) = qt_mesh%cell_vcnr(ctgt,1)
qt_mesh%cell_vcnr(cins,2) = qt_mesh%cell_vemid(ctgt,1)
qt_mesh%cell_vcnr(cins,3) = vtxCC
qt_mesh%cell_vcnr(cins,4) = qt_mesh%cell_vemid(ctgt,4)
qt_mesh%cell_child(ctgt,1) = cins
qt_mesh%cell_level(cins) = qt_mesh%cell_level(ctgt) + 1
cins = cins + 1
if (cins .GE. cm2dopt%Ncell_max) then
    cm2dopt%cm2dfailure = 1
    write(*,'(A)') 'maximum cell count exceeded'
    return 
end if

!child 2
qt_mesh%cell_vcnr(cins,1) = qt_mesh%cell_vemid(ctgt,1)
qt_mesh%cell_vcnr(cins,2) = qt_mesh%cell_vcnr(ctgt,2)
qt_mesh%cell_vcnr(cins,3) = qt_mesh%cell_vemid(ctgt,2)
qt_mesh%cell_vcnr(cins,4) = vtxCC
qt_mesh%cell_child(ctgt,2) = cins
qt_mesh%cell_level(cins) = qt_mesh%cell_level(ctgt) + 1
cins = cins + 1
if (cins .GE. cm2dopt%Ncell_max) then
    cm2dopt%cm2dfailure = 1
    write(*,'(A)') 'maximum cell count exceeded'
    return 
end if

!child 3
qt_mesh%cell_vcnr(cins,1) = vtxCC
qt_mesh%cell_vcnr(cins,2) = qt_mesh%cell_vemid(ctgt,2)
qt_mesh%cell_vcnr(cins,3) = qt_mesh%cell_vcnr(ctgt,3)
qt_mesh%cell_vcnr(cins,4) = qt_mesh%cell_vemid(ctgt,3)
qt_mesh%cell_child(ctgt,3) = cins
qt_mesh%cell_level(cins) = qt_mesh%cell_level(ctgt) + 1
cins = cins + 1
if (cins .GE. cm2dopt%Ncell_max) then
    cm2dopt%cm2dfailure = 1
    write(*,'(A)') 'maximum cell count exceeded'
    return 
end if

!child 4
qt_mesh%cell_vcnr(cins,1) = qt_mesh%cell_vemid(ctgt,4)
qt_mesh%cell_vcnr(cins,2) = vtxCC
qt_mesh%cell_vcnr(cins,3) = qt_mesh%cell_vemid(ctgt,3) 
qt_mesh%cell_vcnr(cins,4) = qt_mesh%cell_vcnr(ctgt,4)
qt_mesh%cell_child(ctgt,4) = cins
qt_mesh%cell_level(cins) = qt_mesh%cell_level(ctgt) + 1
cins = cins + 1
if (cins .GE. cm2dopt%Ncell_max) then
    cm2dopt%cm2dfailure = 1
    write(*,'(A)') 'maximum cell count exceeded'
    return 
end if

!Set and update cell adjacency in this region 
call set_local_adjacency(qt_mesh,ccadj,eopposite,ctgt)

!Flood any edge or face midpoint vertices through the new child cells
do ccf=1,4

    !Child index
    child_idx = qt_mesh%cell_child(ctgt,ccf)

    !Pull adjacent edge midpoint vertices
    do aa=1,4 
        cadj = qt_mesh%cell_adjacent(child_idx,aa)
        if (qt_mesh%cell_vemid(child_idx,aa) == 0) then 
            if (cadj .GT. 0) then 
                if (qt_mesh%cell_level(cadj) == qt_mesh%cell_level(child_idx)) then 
                    qt_mesh%cell_vemid(child_idx,aa) = qt_mesh%cell_vemid(cadj,eopposite(aa))
                end if 
            end if
        end if 
    end do 

    !Push edge midpoint vertices 
    do aa=1,4 
        cadj = qt_mesh%cell_adjacent(child_idx,aa)
        if (cadj .GT. 0) then 
            if (qt_mesh%cell_level(cadj) == qt_mesh%cell_level(child_idx)) then 
                if (qt_mesh%cell_vemid(cadj,eopposite(aa)) == 0) then 
                    qt_mesh%cell_vemid(cadj,eopposite(aa)) = qt_mesh%cell_vemid(child_idx,aa)
                end if
            end if 
        end if 
    end do 
end do 
return 
end subroutine quadtree_cell_refine




!Cell child adjacency subroutine ===========================
subroutine set_local_adjacency(qt_mesh,ccadj,eopposite,cidx)
implicit none

!Variables - Import
integer(in) :: cidx
integer(in) :: eopposite(4),ccadj(4,4)
type(quadtree_data) :: qt_mesh

!Variables - Local
integer(in) :: cc,aa,child_idx,cadj,cadjidx,etgt_op 

!Set child adjacency for each child cell
do cc=1,4
    child_idx = qt_mesh%cell_child(cidx,cc)
    do aa=1,4
        if (ccadj(cc,aa) .GT. 0) then !Parent internal cell
            qt_mesh%cell_adjacent(child_idx,aa) = qt_mesh%cell_child(cidx,ccadj(cc,aa))
        else !Parent external cell
            cadj = qt_mesh%cell_adjacent(cidx,aa) 
            if (cadj .GT. 0) then !Cell
                if (qt_mesh%cell_child(cadj,1) == 0) then !Adjacent cell has no children
                    qt_mesh%cell_adjacent(child_idx,aa) = cadj
                else !Adjacent cell has children

                    !Find index of adjacent child cell adjacent to this new child cell
                    cadjidx = qt_mesh%cell_child(cadj,abs(ccadj(cc,aa)))

                    !Set adjacency of this parents child to the adjacent cells child
                    qt_mesh%cell_adjacent(child_idx,aa) = cadjidx

                    !Opposite face to the targeted child face
                    etgt_op = eopposite(aa)

                    !Set adjacency of adjacent child cadjidx to this parents child on the opposite face
                    qt_mesh%cell_adjacent(cadjidx,etgt_op) = child_idx

                    !Cascade the adjacency of all children of this adjacent cell cadjidx in the direction ftgt to child_idx if they are tagged with adjacency cref on face ftgt currently
                    call cascade_adjacency(qt_mesh,cadjidx,etgt_op,child_idx,cidx)
                end if 
            else !External boundary
                qt_mesh%cell_adjacent(child_idx,aa) = -2
            end if
        end if
    end do 
end do 
end subroutine set_local_adjacency




!Cell child adjacency cascade subroutine ===========================
recursive subroutine cascade_adjacency(qt_mesh,cell2up,etgt,newadj,oldadj)
implicit none 

!Variables - Import
integer(in) :: cell2up,etgt,newadj,oldadj
type(quadtree_data), intent(inout) :: qt_mesh

!Variables - Local
integer(in) :: ch,cell2upN 

!Update current cell cell2up
if (qt_mesh%cell_adjacent(cell2up,etgt) == oldadj) then 
    qt_mesh%cell_adjacent(cell2up,etgt) = newadj
end if 

!Check for and update children of cell cell2up
do ch=1,4
    if (qt_mesh%cell_child(cell2up,ch) .NE. 0) then 
        cell2upN = qt_mesh%cell_child(cell2up,ch)
        call cascade_adjacency(qt_mesh,cell2upN,etgt,newadj,oldadj)
    end if
end do 
return 
end subroutine cascade_adjacency




!Segment - cell overlap testing function ===========================
function segment_cell_ovlp_bool(qt_mesh,surface_mesh,surface_adtree,node_select,nselected,cidx) result(olstat)
implicit none 

!Variables - Import
integer(in) :: olstat,nselected,cidx
integer(in), dimension(:) :: node_select
type(tree_data) :: surface_adtree
type(surface_data) :: surface_mesh
type(quadtree_data) :: qt_mesh

!Variables - Local
integer(in) :: nn,kk,exist 
real(dp) :: vl1(2),vl2(2),v1(2),v2(2),v3(2),v4(2)

!Test
exist = 0 
olstat = 0 
v1(:) = qt_mesh%vtx(qt_mesh%cell_vcnr(cidx,1),:)
v2(:) = qt_mesh%vtx(qt_mesh%cell_vcnr(cidx,2),:)
v3(:) = qt_mesh%vtx(qt_mesh%cell_vcnr(cidx,3),:)
v4(:) = qt_mesh%vtx(qt_mesh%cell_vcnr(cidx,4),:)
do nn=1,nselected
    do kk=1,surface_adtree%tree(node_select(nn))%nentry

        !Verticies of the ends of the surface segment 
        vl1(:) = surface_mesh%vertices(surface_mesh%faces(surface_adtree%tree(node_select(nn))%entry(kk),1),:)
        vl2(:) = surface_mesh%vertices(surface_mesh%faces(surface_adtree%tree(node_select(nn))%entry(kk),2),:)

        !Check for overlap of surface line segment with cell
        exist = seg_rectangle_intersect_bool(v1,v2,v3,v4,vl1,vl2)

        !Exit if overlap found
        if (exist .NE. 0) then 
            olstat = 1
            exit 
        end if
    end do
    if (olstat == 1) then
        exit 
    end if
end do 
return 
end function segment_cell_ovlp_bool




!Segment - cell overlap with surface radius of curvature check testing function ===========================
function segment_cell_ovlp_wrcurv_bool(qt_mesh,surface_mesh,surface_adtree,node_select,nselected,cidx) result(olstat)
implicit none 

!Variables - Import
integer(in) :: olstat,nselected,cidx
integer(in), dimension(:) :: node_select
type(tree_data) :: surface_adtree
type(surface_data) :: surface_mesh
type(quadtree_data) :: qt_mesh

!Variables - Local
integer(in) :: nn,kk,exist,curv_valid 
real(dp) :: sRcurv,cdx,cdy,reflen
real(dp) :: vl1(2),vl2(2),v1(2),v2(2),v3(2),v4(2)
real(dp) :: xmin,xmax,ymin,ymax

!Test
exist = 0 
olstat = 0 
curv_valid = 0
v1(:) = qt_mesh%vtx(qt_mesh%cell_vcnr(cidx,1),:)
v2(:) = qt_mesh%vtx(qt_mesh%cell_vcnr(cidx,2),:)
v3(:) = qt_mesh%vtx(qt_mesh%cell_vcnr(cidx,3),:)
v4(:) = qt_mesh%vtx(qt_mesh%cell_vcnr(cidx,4),:)
cdx = abs(max(v1(1),v2(1),v3(1),v4(1)) - min(v1(1),v2(1),v3(1),v4(1)))
cdy = abs(max(v1(2),v2(2),v3(2),v4(2)) - min(v1(2),v2(2),v3(2),v4(2)))
! reflen = sqrt(cdx*cdy)
reflen = sqrt(cdx*cdx + cdy*cdy)
xmin = min(v1(1),v2(1),v3(1),v4(1))
xmax = max(v1(1),v2(1),v3(1),v4(1))
ymin = min(v1(2),v2(2),v3(2),v4(2))
ymax = max(v1(2),v2(2),v3(2),v4(2))
do nn=1,nselected
    do kk=1,surface_adtree%tree(node_select(nn))%nentry

        !Verticies of the ends of the surface segment 
        vl1(:) = surface_mesh%vertices(surface_mesh%faces(surface_adtree%tree(node_select(nn))%entry(kk),1),:)
        vl2(:) = surface_mesh%vertices(surface_mesh%faces(surface_adtree%tree(node_select(nn))%entry(kk),2),:)

        !Check for overlap of surface line segment with cell
        exist = seg_rectangle_intersect_bool(v1,v2,v3,v4,vl1,vl2)
        
        !Check if surface radius of curvature here is smaller than the cell size
        sRcurv = surface_mesh%face_rcurv(surface_adtree%tree(node_select(nn))%entry(kk))
        if (sRcurv .LE. reflen) then 
            curv_valid = 1
        else
            curv_valid = 0 
        end if 

        ! !Check if either vertex is within the cell and has a surface radius of curvature smaller than the cell size
        ! if ((vl1(1) .GE. xmin) .AND. (vl1(1) .LE. xmax)) then 
        !     if ((vl1(2) .GE. ymin) .AND. (vl1(2) .LE. ymax)) then 
        !         sRcurv = surface_mesh%vtx_rcurv(surface_mesh%faces(surface_adtree%tree(node_select(nn))%entry(kk),1))
        !         if (sRcurv .LE. reflen) then 
        !             curv_valid = 1
        !         end if 
        !     end if 
        ! end if 
        ! if ((vl2(1) .GE. xmin) .AND. (vl2(1) .LE. xmax)) then 
        !     if ((vl2(2) .GE. ymin) .AND. (vl2(2) .LE. ymax)) then 
        !         sRcurv = surface_mesh%vtx_rcurv(surface_mesh%faces(surface_adtree%tree(node_select(nn))%entry(kk),2))
        !         if (sRcurv .LE. reflen) then 
        !             curv_valid = 1
        !         end if 
        !     end if 
        ! end if 

        !Exit if overlap found and curvature condition is met
        if ((exist .NE. 0) .AND. (curv_valid == 1)) then 
            olstat = 1
            exit 
        end if
    end do
    if (olstat == 1) then
        exit 
    end if
end do 
return 
end function segment_cell_ovlp_wrcurv_bool




!Adjacency difference calculation function ===========================
function max_adjacent_dR(qt_mesh,cidx) result(dR)
implicit none 

!Variables - Import
integer(in) :: dR,cidx
type(quadtree_data) :: qt_mesh

!Variables - Local
integer(in) :: kk,lc 
integer(in) :: deltaR(4)

!Check adjacent cells 
lc = qt_mesh%cell_level(cidx)
do kk=1,4
    if (qt_mesh%cell_adjacent(cidx,kk) .GT. 0) then
        deltaR(kk) = qt_mesh%cell_level(qt_mesh%cell_adjacent(cidx,kk)) - lc
    else
        deltaR(kk) = 0 
    end if
end do 
dR = maxval(deltaR(:))
return 
end function max_adjacent_dR




!Quadtreee trim to geometry function ===========================
subroutine quadtree_mesh_trim2geom(cell_keep,vtx_external,vtx_type1,qt_mesh,surface_mesh,surface_adtree,cm2dopt) 
implicit none 

!System data
type(quadtree_data) :: qt_mesh
type(surface_data) :: surface_mesh
type(cm2d_options) :: cm2dopt
type(tree_data) :: surface_adtree

!Variables - Import
integer(in), dimension(:), allocatable :: cell_keep,vtx_external,vtx_type1

!Variables - Local 
integer(in) :: cc,ee,nn,kk,ii,aa,vv
integer(in) :: Npert,nselected,Noverlap,v1,v2,segval,edgenew,intnew,segselect,edgetgt,stransstat,cadj
integer(in) :: Nfront,NfrontN,Nexternal,Ninternal,cellC,cellA,etgt_op,Nflooditer,ctgt_adj,etgt_adj,Nsubcedge_vtx
integer(in) :: eopposite(4),cell_nc(qt_mesh%cins-1),node_select(surface_adtree%nnode),subcedge_vtx(cm2dopt%NintEmax)
integer(in) :: frontC(qt_mesh%cins-1),frontN(qt_mesh%cins-1),cell_external_base(qt_mesh%cins-1)
integer(in) :: celledges(4,2),edge_visit(qt_mesh%cins-1,4)
integer(in) :: seg_intersect_edge(cm2dopt%NintEmax),seg_intersect_edge_odr(cm2dopt%NintEmax),int_selected(cm2dopt%NintEmax)
real(dp) :: transpert,cpadSZ,zxmin,zxmax,zymin,zymax,zzmin,zzmax,segnorm
real(dp) :: ve1(2),ve2(2),segdir(2),vl1(2),vl2(2),vint(2),segnormal(2),segdist(cm2dopt%NintEmax)
real(dp) :: seg_intersect(cm2dopt%NintEmax,2),seg_intersect_odr(cm2dopt%NintEmax,2)

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

!Define intersection transition perturbation 
transpert = cm2dopt%intcointol 

!Allocate cell_keep, vtx_external and vtx_type1
allocate(cell_keep(qt_mesh%cins-1))
allocate(vtx_external(qt_mesh%vins-1))
allocate(vtx_type1(qt_mesh%vins-1))

!Initialise cell keep states 0 = removed || 1 = geometry overlap || 2 = geometry external || 3 = geometry internal 
cell_keep(:) = 0 

!Identify cells with no children 
cell_nc(:) = 0 
do cc=1,qt_mesh%cins-1
    if (qt_mesh%cell_child(cc,1) == 0) then !Has no children
        cell_nc(cc) = 1
    end if 
end do

!Display 
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '--> tie-breaking vertices that lie on geometry surfaces'
end if

!Perturb vertices that lie exactly on geometry surfaces 
do cc=1,5
    call perturb_vertices_on_surfaces(Npert,qt_mesh,surface_mesh,cell_nc,surface_adtree,cm2dopt) 
    !print *, Npert
    if (Npert == 0) then 
        exit 
    end if 
end do 

!Display 
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '--> identifying cell states with respect to the surface geometry'
end if

!Identify all valid cells with no children that overlap the geometry boundary and clasify vertices on these overlap cells as geometry internal or external 
zxmin = 0.0d0 
zxmax = 0.0d0 
zymin = 0.0d0 
zymax = 0.0d0 
zzmin = 0.0d0 
zzmax = 0.0d0
Noverlap = 0 
nselected = 0 
node_select(:) = 0
vtx_external(:) = 0
edge_visit(:,:) = 0 
seg_intersect(:,:) = 0.0d0  
seg_intersect_edge(:) = 0.0d0  
do cc=1,qt_mesh%cins-1
    if (cell_nc(cc) == 1) then 
        
        !Padding size 
        cpadSZ = 2.0d0*cm2dopt%far_field_bound/(2.0d0**(qt_mesh%cell_level(cc) - 1))

        !Intersection bounding box
        zxmin = qt_mesh%vtx(qt_mesh%cell_vcnr(cc,1),1) - cpadSZ*cm2dopt%ADTpadding !tgt bounding box -> xmin
        zxmax = qt_mesh%vtx(qt_mesh%cell_vcnr(cc,3),1) + cpadSZ*cm2dopt%ADTpadding !tgt bounding box -> xmax
        zymin = qt_mesh%vtx(qt_mesh%cell_vcnr(cc,1),2) - cpadSZ*cm2dopt%ADTpadding !tgt bounding box -> ymin
        zymax = qt_mesh%vtx(qt_mesh%cell_vcnr(cc,3),2) + cpadSZ*cm2dopt%ADTpadding !tgt bounding box -> ymax

        !Identify any segment bounding boxes that may overlap the cell 
        call search_ADtree(nselected,node_select,surface_adtree,zxmin,zxmax,zymin,zymax,zzmin,zzmax)

        !Check for and tag geometry surface overlap 
        if (segment_cell_ovlp_bool(qt_mesh,surface_mesh,surface_adtree,node_select,nselected,cc) == 1) then 
            cell_keep(cc) = 1
        end if

        !Classify vertices on this cell
        if (cell_keep(cc) == 1) then 

            !Search edges
            do ee=1,4
                if (edge_visit(cc,ee) == 0) then 

                    !Set edge to visited
                    edge_visit(cc,ee) = 1 

                    !Edge vertices
                    v1 = qt_mesh%cell_vcnr(cc,celledges(ee,1))
                    v2 = qt_mesh%cell_vcnr(cc,celledges(ee,2))
                    ve1(:) = qt_mesh%vtx(v1,:)
                    ve2(:) = qt_mesh%vtx(v2,:)

                    !Edge direction
                    segdir(:) = ve2(:) - ve1(:)
                    segnorm = norm2(segdir(:))

                    !Check for intersections on this edge 
                    Noverlap = 0 
                    seg_intersect_edge(:) = 0 
                    do nn=1,nselected
                        do kk=1,surface_adtree%tree(node_select(nn))%nentry

                            !Verticies of the ends of the surface segment 
                            vl1(:) = surface_mesh%vertices(surface_mesh%faces(surface_adtree%tree(node_select(nn))%entry(kk),1),:)
                            vl2(:) = surface_mesh%vertices(surface_mesh%faces(surface_adtree%tree(node_select(nn))%entry(kk),2),:)

                            !If potential intersection
                            segval = seg_seg_intersect_bool(vl1,vl2,ve1,ve2)
                            if (segval .NE. 0) then !Only select where intersection is properly defined and there is a defined crossing point (dont need to find all cases here but must be correct)

                                !Find intersection location (is within edge and segment here)
                                vint = line_line_intersection_loc_inl1(ve1,ve2,vl1,vl2)
              
                                !Test against existing intersections to ensure it is new and valid --- 
                                !Check if a new surface segment
                                edgenew = 1 
                                do ii=1,Noverlap
                                    if (surface_adtree%tree(node_select(nn))%entry(kk) == seg_intersect_edge(Noverlap)) then 
                                        edgenew = 0 
                                    end if
                                end do  
        
                                !Check if co-incidint with another intersection on this edge 
                                intnew = 1
                                do ii=1,Noverlap
                                    if (norm2(vint(:)-seg_intersect(ii,:)) .LE. cm2dopt%intcointol*segnorm) then
                                        intnew = 0 
                                    end if 
                                end do  

                                !Add intersection if valid and new 
                                if ((edgenew == 1) .AND. (intnew == 1)) then 
                                    Noverlap = Noverlap + 1
                                    if (Noverlap .GT. cm2dopt%NintEmax) then 
                                        cm2dopt%cm2dfailure = 1
                                        print *, '** maximum number of edge-geometry intersections exceeded on edge ',ee
                                        return 
                                    end if 
                                    seg_intersect(Noverlap,:) = vint(:)
                                    seg_intersect_edge(Noverlap) = surface_adtree%tree(node_select(nn))%entry(kk)
                                end if
                            end if 
                        end do 
                    end do 

                    !Order from start to end of edge if required
                    if (Noverlap .GE. 2) then 
                        segdist(:) = 0 
                        do nn=1,Noverlap
                            segdist(nn) = norm2(seg_intersect(nn,:) - ve1(:))
                        end do 
                        int_selected(:) = 0 
                        do nn=1,Noverlap

                            !Select intersection at maximum distance
                            segselect = maxloc(segdist(1:Noverlap),1,mask = int_selected(1:Noverlap) == 0)

                            !Store 
                            seg_intersect_edge_odr(Noverlap-nn+1) = seg_intersect_edge(segselect)
                            seg_intersect_odr(Noverlap-nn+1,:) = seg_intersect(segselect,:)

                            !Set distance to zero and set to visited
                            int_selected(segselect) = 1
                            segdist(segselect) = 0.0d0 
                        end do 
                    else
                        seg_intersect_edge_odr(:) = seg_intersect_edge(:)
                        seg_intersect_odr(:,:) = seg_intersect(:,:)
                    end if 

                    !If any intersections, then classify each end of this edge     
                    if (Noverlap .GT. 0) then 

                        !Classify end 1
                        edgetgt = seg_intersect_edge_odr(1)
                        vl1(:) = surface_mesh%vertices(surface_mesh%faces(edgetgt,1),:)
                        vl2(:) = surface_mesh%vertices(surface_mesh%faces(edgetgt,2),:)
                        segnormal(1) = vl1(2) - vl2(2)
                        segnormal(2) = vl2(1) - vl1(1)
                        if (dot_product(segnormal,segdir) .GT. 0.0d0) then !Exit
                            stransstat = -1
                        elseif (dot_product(segnormal,segdir) .LT. 0.0d0) then !Entry
                            stransstat = 1
                        else !No transition 
                            stransstat = 0
                        end if 
                        if (stransstat == 1) then !Entry so v1 external
                            vtx_external(v1) = 1
                        end if

                        !Classify end 2
                        edgetgt = seg_intersect_edge_odr(Noverlap)
                        vl1(:) = surface_mesh%vertices(surface_mesh%faces(edgetgt,1),:)
                        vl2(:) = surface_mesh%vertices(surface_mesh%faces(edgetgt,2),:)
                        segnormal(1) = vl1(2) - vl2(2)
                        segnormal(2) = vl2(1) - vl1(1)
                        if (dot_product(segnormal,segdir) .GT. 0.0d0) then !Exit
                            stransstat = -1
                        elseif (dot_product(segnormal,segdir) .LT. 0.0d0) then !Entry
                            stransstat = 1
                        else !No transition 
                            stransstat = 0
                        end if 
                        if (stransstat == -1) then !Exit so v2 external
                            vtx_external(v2) = 1
                        end if
                    end if

                    !Flood visited edge status
                    do aa=1,4
                        cadj = qt_mesh%cell_adjacent(cc,aa)
                        if (cadj .GT. 0) then 
                            if (cell_nc(cadj) == 1) then 
                                if (qt_mesh%cell_level(cadj) == qt_mesh%cell_level(cc)) then 
                                    edge_visit(cadj,eopposite(aa)) = edge_visit(cc,eopposite(aa))
                                end if 
                            end if 
                        end if 
                    end do
                end if
            end do 
        end if 
    end if 
end do 

!Tag as external base cells any with a vertex tagged as external that are not a cell_keep = 1 cell
cell_external_base(:) = 0 
do cc=1,qt_mesh%cins-1
    if ((cell_nc(cc) == 1) .AND. (cell_keep(cc) == 0)) then 
        if (maxval(vtx_external(qt_mesh%cell_vcnr(cc,:))) == 1) then 
            cell_external_base(cc) = 1
        end if
    end if
end do 

!Initialise flood of geometry external cells 
Nfront = 0
Nexternal = 0
frontN(:) = 0
frontC(:) = 0
do cc=1,qt_mesh%cins-1
    if (cell_external_base(cc) == 1) then 
        cell_keep(cc) = 2
        Nfront = Nfront + 1
        Nexternal = Nexternal + 1
        frontC(Nfront) = cc 
    end if
end do 

!Flood external cells 
Nflooditer = 0 
do nn=1,qt_mesh%cins-1 !Flood iterations 

    !Reset new front count 
    NfrontN = 0

    !Flood to adjacent cells on the front 
    do ii=1,Nfront

        !Current cell
        cellC = frontC(ii)

        !Each adjacent cell of this cell 
        do ee=1,4
            cellA = qt_mesh%cell_adjacent(cellC,ee)
            if (cellA .GT. 0) then 
                if (cell_nc(cellA) == 1) then !If valid cell 
                    if (cell_keep(cellA) == 0) then !Not yet visited -> tag as external and add to the front 
                        NfrontN = NfrontN + 1
                        cell_keep(cellA) = 2
                        Nexternal = Nexternal + 1
                        frontN(NfrontN) = cellA
                    end if
                else !Has children -> check and tag valid children adjacent to this edge 
                    etgt_op = eopposite(ee)
                    call find_adjacent_children(cell_keep,frontN,NfrontN,Nexternal,cellC,cellA,etgt_op,qt_mesh,cell_nc)
                end if
            end if
        end do
    end do
    
    !Exit if no new cells found
    if (NfrontN == 0) then 
        Nflooditer = nn
        exit 
    end if

    !Update flood items for next iteration
    Nfront = NfrontN 
    frontC(1:NfrontN) = frontN(1:NfrontN)
end do 

!Set all vertices in external cells to external 
do cc=1,qt_mesh%cins-1
    if (cell_keep(cc) == 2) then 

        !Set base cell vertices 
        vtx_external(qt_mesh%cell_vcnr(cc,:)) = 1
        do ee=1,4
            if (qt_mesh%cell_vemid(cc,ee) .NE. 0) then 
                vtx_external(qt_mesh%cell_vemid(cc,ee)) = 1
            end if 
        end do 

        !Set any edge mid vertices from adjacent cells as external on this edge 
        Nsubcedge_vtx = 0 
        subcedge_vtx(:) = 0 
        do ee=1,4 
            ctgt_adj = qt_mesh%cell_adjacent(cc,ee)
            etgt_adj = eopposite(ee)
            if (ctgt_adj .GT. 0) then 
                if (qt_mesh%cell_level(ctgt_adj) == qt_mesh%cell_level(cc)) then 
                    call get_qtmesh_edge_mid_vertices(subcedge_vtx,Nsubcedge_vtx,ctgt_adj,etgt_adj,qt_mesh,celledges)
                end if 
            end if 
        end do 
        if (Nsubcedge_vtx .NE. 0) then 
            do vv=1,Nsubcedge_vtx
                vtx_external(subcedge_vtx(vv)) = 1
            end do 
        end if
    end if 
end do 

!Tag all vertices that are part of a type 1 intersecting cell
vtx_type1(:) = 0 
do cc=1,qt_mesh%cins-1
    if (cell_nc(cc) == 1) then 
        if (cell_keep(cc) == 1) then 
            vtx_type1(qt_mesh%cell_vcnr(cc,:)) = 1
        end if 
    end if 
end do 

!Display
if (cm2dopt%dispt == 1) then
    write(*,'(A,I0,A,I0,A)') '    {',Nexternal,' external cells identified after ',Nflooditer,' iterations}'
end if

!Tag all geometry internal cells 
Ninternal = 0 
do cc=1,qt_mesh%cins-1
    if (cell_nc(cc) == 1) then 
        if (cell_keep(cc) == 0) then 
            Ninternal = Ninternal + 1
            cell_keep(cc) = 3
        end if
    end if
end do 

!Display
if (cm2dopt%dispt == 1) then
    write(*,'(A,I0,A)') '    {',Ninternal,' internal cells identified}'
end if

!Debug write cells of a type -------
! open(11,file='io/celloftype')
! do cc=1,qt_mesh%cins-1
!     if (cell_nc(cc) == 1) then 
!         if (cell_keep(cc) == 2) then  
!         ! if (cell_external_base(cc) == 1) then 
!             write(11,*) qt_mesh%vtx(qt_mesh%cell_vcnr(cc,1),:)
!             write(11,*) qt_mesh%vtx(qt_mesh%cell_vcnr(cc,2),:)
!             write(11,*) qt_mesh%vtx(qt_mesh%cell_vcnr(cc,3),:)
!             write(11,*) qt_mesh%vtx(qt_mesh%cell_vcnr(cc,4),:)
!         end if
!     end if
! end do 
! close(11)

!Debug write reference external vertices -------
! open(11,file='io/vtxexternal')
! do cc=1,qt_mesh%vins-1
!     if (vtx_external(cc) == 1) then 
!         write(11,*) qt_mesh%vtx(cc,:)
!     end if
! end do 
! close(11)
return 
end subroutine quadtree_mesh_trim2geom




!Surface cooincident vertex perturbation subroutine ===========================
subroutine perturb_vertices_on_surfaces(Npert,qt_mesh,surface_mesh,cell_nc,surface_adtree,cm2dopt) 
use ieee_arithmetic
implicit none 

!System data
type(quadtree_data) :: qt_mesh
type(surface_data) :: surface_mesh
type(cm2d_options) :: cm2dopt
type(tree_data) :: surface_adtree

!Variables - Import
integer(in) :: Npert
integer(in) :: cell_nc(qt_mesh%cins-1)

!Variables - Local 
integer(in) :: cc,vv,nn,kk,nselected,vincident 
integer(in) :: node_select(surface_adtree%nnode),vtx_updated(qt_mesh%vins-1)
real(dp) :: cpadSZ,surf_tol,zxmin,zxmax,zymin,zymax,zzmin,zzmax
real(dp) :: vtxb(2),vl1(2),vl2(2),segnorm(2),vid(3)

!Set proximity tollerance as fraction of cell size
surf_tol = 2.0d0*cm2dopt%intcointol 

!Check and perturb vertices 
zxmin = 0.0d0 
zxmax = 0.0d0 
zymin = 0.0d0 
zymax = 0.0d0 
zzmin = 0.0d0 
zzmax = 0.0d0
nselected = 0  
node_select(:) = 0 
vtx_updated(:) = 0 
!open(11,file=cm2dopt%iopath//'vtx2pert.dat')
do cc=1,qt_mesh%cins-1
    if (cell_nc(cc) == 1) then 
        
        !Padding size 
        cpadSZ = 2.0d0*cm2dopt%far_field_bound/(2.0d0**(qt_mesh%cell_level(cc) - 1))

        !Intersection bounding box
        zxmin = qt_mesh%vtx(qt_mesh%cell_vcnr(cc,1),1) - cpadSZ*cm2dopt%ADTpadding !tgt bounding box -> xmin
        zxmax = qt_mesh%vtx(qt_mesh%cell_vcnr(cc,3),1) + cpadSZ*cm2dopt%ADTpadding !tgt bounding box -> xmax
        zymin = qt_mesh%vtx(qt_mesh%cell_vcnr(cc,1),2) - cpadSZ*cm2dopt%ADTpadding !tgt bounding box -> ymin
        zymax = qt_mesh%vtx(qt_mesh%cell_vcnr(cc,3),2) + cpadSZ*cm2dopt%ADTpadding !tgt bounding box -> ymax

        !Identify any segment bounding boxes that may overlap the cell 
        call search_ADtree(nselected,node_select,surface_adtree,zxmin,zxmax,zymin,zymax,zzmin,zzmax)

        !Check each vertex on this cell for proximity 
        do vv=1,4
            if (vtx_updated(qt_mesh%cell_vcnr(cc,vv)) == 0) then 
                vincident = 0 
                if (nselected .GT. 0) then 
                    do nn=1,nselected
                        do kk=1,surface_adtree%tree(node_select(nn))%nentry

                            !Vertex location
                            vtxb(:) = qt_mesh%vtx(qt_mesh%cell_vcnr(cc,vv),:)

                            !Verticies on the ends of the surface segment 
                            vl1(:) = surface_mesh%vertices(surface_mesh%faces(surface_adtree%tree(node_select(nn))%entry(kk),1),:)
                            vl2(:) = surface_mesh%vertices(surface_mesh%faces(surface_adtree%tree(node_select(nn))%entry(kk),2),:)

                            !Check if vertex is close to this segment 
                            if (norm2(vl2(:) - vl1(:)) .NE. 0.0d0) then 
                                vid = min_dist_point_to_edge(vl1,vl2,vtxb)
                            else
                                vid(3) = ieee_value(1.0d0,IEEE_POSITIVE_INF)
                            end if 
            
                            !If vertex is within tollerance then perturb
                            ! if (vid(3) .LE. 100.0d0*cm2dopt%intcointol) then 
                            if (vid(3) .LE. surf_tol*cpadSZ) then 
                                vincident = 1
                                segnorm(1) = vl1(2) - vl2(2)
                                segnorm(2) = vl2(1) - vl1(1)
                                exit 
                            end if 
                        end do
                        if (vincident == 1) then 
                            exit 
                        end if 
                    end do 
                end if 

                !If vetex is on surface to tollerance 
                if (vincident == 1) then 
                    if (cm2dopt%meshinout == 'out') then !External mesh -> move inside geometry
                        if (norm2(segnorm(:)) .NE. 0.0d0) then 
                            qt_mesh%vtx(qt_mesh%cell_vcnr(cc,vv),:) = qt_mesh%vtx(qt_mesh%cell_vcnr(cc,vv),:) - &
                            (segnorm(:)/norm2(segnorm(:)))*cpadSZ*surf_tol
                        else

                        end if 
                    elseif (cm2dopt%meshinout == 'in') then !Internal mesh -> move outside geometry
                        if (norm2(segnorm(:)) .NE. 0.0d0) then 
                            qt_mesh%vtx(qt_mesh%cell_vcnr(cc,vv),:) = qt_mesh%vtx(qt_mesh%cell_vcnr(cc,vv),:) + &
                            (segnorm(:)/norm2(segnorm(:)))*cpadSZ*surf_tol
                        else

                        end if 
                    end if 
                    vtx_updated(qt_mesh%cell_vcnr(cc,vv)) = 1
                    ! write(11,*) vtxb(:) 
                    ! write(11,*) qt_mesh%vtx(qt_mesh%cell_vcnr(cc,vv),:)
                end if 
            end if 
        end do 
    end if 
end do 
!close(11)

!Set number of updated vertices 
Npert = sum(vtx_updated(:))
return 
end subroutine perturb_vertices_on_surfaces




!Identify adjacent child cells cascade subroutine ===========================
recursive subroutine find_adjacent_children(cell_keep,frontN,NfrontN,Nexternal,cell_tgt,cell2check,etgt_op,qt_mesh,cell_nc)
implicit none 

!System data
type(quadtree_data), intent(in) :: qt_mesh

!Variables - Import
integer(in), intent(inout) :: NfrontN,Nexternal
integer(in) :: cell_tgt,cell2check,etgt_op
integer(in) :: cell_nc(qt_mesh%cins-1)
integer(in), intent(inout) :: cell_keep(qt_mesh%cins-1),frontN(qt_mesh%cins-1)

!Variables - Local 
integer(in) :: ch,child_idx

!Check children of cell cell2check for those that are adjacent to cell_tgt that have no children 
do ch=1,4
    if (qt_mesh%cell_child(cell2check,ch) .NE. 0) then 
        child_idx = qt_mesh%cell_child(cell2check,ch)
        if (qt_mesh%cell_adjacent(child_idx,etgt_op) == cell_tgt) then !If adjacent to current target cell 
            if (cell_nc(child_idx) == 1) then !If cell child_idx has no children 
                if (cell_keep(child_idx) == 0) then !Not yet visited -> hence add to new front and set as external
                    NfrontN = NfrontN + 1
                    Nexternal = Nexternal + 1
                    cell_keep(child_idx) = 2
                    frontN(NfrontN) = child_idx
                end if
            else !Search this cells children for valid cells if any -> recursive call 
                call find_adjacent_children(cell_keep,frontN,NfrontN,Nexternal,cell_tgt,child_idx,etgt_op,qt_mesh,cell_nc)
            end if
        end if
    end if 
end do 
return 
end subroutine find_adjacent_children




!Subedge vertex identification subroutine ===========================
recursive subroutine get_qtmesh_edge_mid_vertices(subcedge_vtx,Nsubcedge_vtx,ctgt,etgt,qt_mesh,edges)
implicit none 

!Variables - Import 
integer(in), intent(inout) :: ctgt,etgt,Nsubcedge_vtx
integer(in) :: edges(4,2)
integer(in), dimension(:) :: subcedge_vtx
type(quadtree_data) :: qt_mesh

!Variables - Local 
integer(in) :: cc,vv,vtgt,exist,c_child,etgt_v1,etgt_v2,echld_v1,echld_v2  

!Check if target edge contains a midpoint vertex
if (qt_mesh%cell_vemid(ctgt,etgt) .NE. 0) then 

    !Target midpoint vertex
    vtgt = qt_mesh%cell_vemid(ctgt,etgt)

    !Check if vertex is new 
    exist = 0 
    do vv=1,Nsubcedge_vtx
        if (subcedge_vtx(vv) == vtgt) then 
            exist = 1 
            exit 
        end if 
    end do 

    !Add this edge midpoint vertex to the list if new
    if (exist == 0) then 
        Nsubcedge_vtx = Nsubcedge_vtx + 1
        subcedge_vtx(Nsubcedge_vtx) = vtgt
    end if 

    !If this cell has children then check the target edge in the children that contain either vertex of etgt
    if (qt_mesh%cell_child(ctgt,1) .NE. 0) then 

        !Vertices in the current cell on the target edge 
        etgt_v1 = qt_mesh%cell_vcnr(ctgt,edges(etgt,1))
        etgt_v2 = qt_mesh%cell_vcnr(ctgt,edges(etgt,2)) 
        
        !Check children
        do cc=1,4

            !Child cell 
            c_child = qt_mesh%cell_child(ctgt,cc)

            !Target edge ends in child cell
            echld_v1 = qt_mesh%cell_vcnr(c_child,edges(etgt,1))
            echld_v2 = qt_mesh%cell_vcnr(c_child,edges(etgt,2))

            !Check this cell and edge if either vertex is shared with the parent on the targeted edge
            if ((echld_v1 == etgt_v1) .OR. (echld_v2 == etgt_v2) .OR. (echld_v1 == etgt_v2) .OR. (echld_v2 == etgt_v1)) then 
                call get_qtmesh_edge_mid_vertices(subcedge_vtx,Nsubcedge_vtx,c_child,etgt,qt_mesh,edges)
            end if 
        end do 
    end if 
end if 
return 
end subroutine get_qtmesh_edge_mid_vertices


end module cellmesh2d_quadtree_mod