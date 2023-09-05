!Module Providing Routines For AD_Tree Construction And Searching For Candidate Line-Object Intersections 
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 4.2
!Updated 29-06-2023

!Module 
module cellmesh2d_adtree_mod
use cellmesh2d_data_mod
implicit none

!ADTree derived types =>
!Tree structure
type ADtree_data
    integer(in) :: nentry,cdntindex
    integer(in) :: children(2)
    real(dp) :: median
    integer(in), dimension(:), allocatable :: entry
    real(dp), dimension(:), allocatable :: ai,bi
end type ADtree_data

!Tree container 
type tree_data 
    integer(in) :: nnode,ndim,max_depth
    real(dp), dimension(:,:), allocatable :: tgt_ab
    type(ADtree_data), dimension(:), allocatable :: tree
end type tree_data 

!Subroutines =>
contains


!Interface to construct tree subroutine =================================================
subroutine build_ADtree(tree_current,ndim,max_depth,node_minDIVsize,tvtx,dpad,disp_toggle)
implicit none 

!Variables - Import
integer(in) :: ndim,max_depth,node_minDIVsize,disp_toggle
real(dp) :: dpad
real(dp), dimension(:,:) :: tvtx
type(tree_data) :: tree_current

!Variables - Local 

!Assign tree data
tree_current%ndim = ndim
tree_current%max_depth = max_depth

!Construct tree 
call ad_tree_construct(tree_current%tree,tree_current%nnode,ndim,max_depth,node_minDIVsize,tvtx,disp_toggle)

!Set search region target box 
if (ndim == 4) then 

    !Allocate
    allocate(tree_current%tgt_ab(4,2))

    !Construct base format for target region bounding box to search adtree
    !ai
    tree_current%tgt_ab(1,1) = minval(tvtx(:,1))
    tree_current%tgt_ab(2,1) = minval(tvtx(:,2))
    tree_current%tgt_ab(3,1) = 0.0d0 !tgt bounding box -> xmin
    tree_current%tgt_ab(4,1) = 0.0d0 !tgt bounding box -> ymin

    !bi
    tree_current%tgt_ab(1,2) = 0.0d0 !tgt bounding box -> xmax
    tree_current%tgt_ab(2,2) = 0.0d0 !tgt bounding box -> ymax
    tree_current%tgt_ab(3,2) = maxval(tvtx(:,3)) 
    tree_current%tgt_ab(4,2) = maxval(tvtx(:,4)) 
elseif (ndim == 6) then 

    !Allocate
    allocate(tree_current%tgt_ab(6,2))

    !Construct base format for target region bounding box
    !ai
    tree_current%tgt_ab(1,1) = minval(tvtx(:,1)) - dpad
    tree_current%tgt_ab(2,1) = minval(tvtx(:,2)) - dpad
    tree_current%tgt_ab(3,1) = minval(tvtx(:,3)) - dpad
    tree_current%tgt_ab(4,1) = 0.0d0 !tgt bounding box -> xmin
    tree_current%tgt_ab(5,1) = 0.0d0 !tgt bounding box -> ymin
    tree_current%tgt_ab(6,1) = 0.0d0 !tgt bounding box -> zmin

    !bi
    tree_current%tgt_ab(1,2) = 0.0d0 !tgt bounding box -> xmax
    tree_current%tgt_ab(2,2) = 0.0d0 !tgt bounding box -> ymax
    tree_current%tgt_ab(3,2) = 0.0d0 !tgt bounding box -> zmax
    tree_current%tgt_ab(4,2) = maxval(tvtx(:,4)) + dpad
    tree_current%tgt_ab(5,2) = maxval(tvtx(:,5)) + dpad
    tree_current%tgt_ab(6,2) = maxval(tvtx(:,6)) + dpad
end if
return 
end subroutine build_ADtree




!Interface to search tree subroutine =================================================
subroutine search_ADtree(nselected,node_select,tree_current,zxmin,zxmax,zymin,zymax,zzmin,zzmax)
implicit none 

!Tree
type(tree_data) :: tree_current

!Variables - Import
integer(in) :: nselected
integer(in) :: node_select(tree_current%nnode)
real(dp) :: zxmin,zxmax,zymin,zymax,zzmin,zzmax

!Dimension switch 
if (tree_current%ndim == 4) then 

    !Set target region bounds 
    tree_current%tgt_ab(3,1) = zxmin !tgt bounding box -> xmin
    tree_current%tgt_ab(4,1) = zymin !tgt bounding box -> ymin
    tree_current%tgt_ab(1,2) = zxmax !tgt bounding box -> xmax
    tree_current%tgt_ab(2,2) = zymax !tgt bounding box -> ymax

    !Search tree for nodes overlapping this region 
    call ad_tree_search4(node_select,nselected,tree_current%tree,tree_current%ndim,tree_current%max_depth,&
                         tree_current%tgt_ab,tree_current%nnode)
elseif (tree_current%ndim == 6) then 

    !Set target region bounds 
    tree_current%tgt_ab(4,1) = zxmin !tgt bounding box -> xmin
    tree_current%tgt_ab(5,1) = zymin !tgt bounding box -> ymin
    tree_current%tgt_ab(6,1) = zzmin !tgt bounding box -> zmin
    tree_current%tgt_ab(1,2) = zxmax !tgt bounding box -> xmax
    tree_current%tgt_ab(2,2) = zymax !tgt bounding box -> ymax
    tree_current%tgt_ab(3,2) = zzmax !tgt bounding box -> zmax

    !Search tree for nodes overlapping this region 
    call ad_tree_search6(node_select,nselected,tree_current%tree,tree_current%ndim,tree_current%max_depth,&
                         tree_current%tgt_ab,tree_current%nnode)
end if 
return 
end subroutine search_ADtree




!Interface to search tree for closest nodes to a point subroutine =================================================
subroutine search_ADtree_closest_nodes2point(nselected,node_select,tree_current,vbase,distref,lpad)
implicit none 

!Tree
type(tree_data) :: tree_current
    
!Variables - Import
integer(in) :: nselected
integer(in) :: node_select(tree_current%nnode)
real(dp) :: distref,lpad
real(dp), dimension(:) :: vbase

!Search tree for the nodes closest to vbase
if (tree_current%ndim == 4) then 

    !Set target vertex padded bounds 
    tree_current%tgt_ab(3,1) = vbase(1) - lpad !tgt bounding box -> xmin
    tree_current%tgt_ab(4,1) = vbase(2) - lpad !tgt bounding box -> ymin
    tree_current%tgt_ab(1,2) = vbase(1) + lpad !tgt bounding box -> xmax
    tree_current%tgt_ab(2,2) = vbase(2) + lpad !tgt bounding box -> ymax

    !Search tree
    call ad_tree_search_closest_nodes_2point4(vbase,distref,node_select,nselected,tree_current%tree,tree_current%ndim,&
                                              tree_current%max_depth,tree_current%tgt_ab,tree_current%nnode)
elseif (tree_current%ndim == 6) then 

end if 
return 
end subroutine search_ADtree_closest_nodes2point




!Interface to search tree for closest nodes to a zone subroutine =================================================
subroutine search_ADtree_closest_nodes2zone(nselected,node_select,tree_current,vz1,vz2,vz3,vz4,distref)
implicit none 

!Tree
type(tree_data) :: tree_current
    
!Variables - Import
integer(in) :: nselected
integer(in) :: node_select(tree_current%nnode)
real(dp) :: distref
real(dp), dimension(:) :: vz1,vz2,vz3,vz4

!Search tree for the nodes closest to vbase
if (tree_current%ndim == 4) then 
    call ad_tree_search_closest_nodes_2zone4(vz1,vz2,vz3,vz4,distref,node_select,nselected,tree_current%tree,tree_current%ndim,&
                                             tree_current%max_depth,tree_current%nnode)
elseif (tree_current%ndim == 6) then 

end if 
return 
end subroutine search_ADtree_closest_nodes2zone




!Interface to search tree for nodes overlapping an arbitrary rectangle subroutine =================================================
subroutine search_ADtree_rec_ovlp(nselected,node_select,tree_current,clv1,clv2,recwidth)
implicit none 

!Tree
type(tree_data) :: tree_current
    
!Variables - Import
integer(in) :: nselected
integer(in) :: node_select(tree_current%nnode)
real(dp) :: recwidth
real(dp), dimension(:) :: clv1,clv2

!Variables - Local 
real(dp), dimension(:), allocatable :: rv1,rv2,rv3,rv4

!Search tree for the nodes closest to vbase
if (tree_current%ndim == 4) then 

    !Build rectangle from centreline vertices   
    call rec_from_centrline_2d(rv1,rv2,rv3,rv4,clv1,clv2,recwidth)

    !Search for overlaps 
    call ad_tree_search_rectangle_ovlp4(node_select,nselected,rv1,rv2,rv3,rv4,tree_current%tree,tree_current%ndim,&
                                        tree_current%max_depth,tree_current%nnode)
elseif (tree_current%ndim == 6) then 

end if 
return 
end subroutine search_ADtree_rec_ovlp




!Subroutine To Construct an AD Tree In ndim Dimensions =================================================
subroutine ad_tree_construct(adtreeI,Nnode,ndim,max_depth,node_minDIVsize,tvtx,disp_toggle)
implicit none 

!Variables - Import
integer(in) :: ndim,max_depth,Nnode,node_minDIVsize,disp_toggle
real(dp), dimension(:,:) :: tvtx
type(ADtree_data), dimension(:), allocatable :: adtreeI

!Variables - Local
integer(in) :: ii,rr,cc,depth,cindex,nvtx,Nchild1,Nchild2,childins,valins
integer(in) :: n_child_node_C,n_child_node_N
integer(in), dimension(:), allocatable :: child_sort,current_node,new_node
real(dp) :: parent_median
type(ADtree_data), dimension(:), allocatable :: adtree

!Properties
nvtx = size(tvtx,1)
allocate(child_sort(nvtx))
allocate(current_node(nvtx))
allocate(new_node(nvtx))

!Allocate tree array
allocate(adtree(2*nvtx))

!Initialise root of tree
current_node(:) = 0
current_node(1) = 1
allocate(adtree(1)%entry(nvtx))
allocate(adtree(1)%ai(ndim))
allocate(adtree(1)%bi(ndim))
do ii=1,nvtx
    adtree(1)%entry(ii) = ii
end do
adtree(1)%nentry = nvtx

!Construct global bounding box
do ii=1,ndim
    adtree(1)%ai(ii) = minval(tvtx(:,ii)) !Minimum of ndim axes
    adtree(1)%bi(ii) = maxval(tvtx(:,ii)) !Maximum of ndim axes
end do

!Construct tree
depth = 0
n_child_node_C = 1
childins = 2
do rr=1,max_depth*ndim

    !Set target coordinate index
    cindex = 1 + mod(depth,ndim)

    !Divide each node along its target coordinate cindex
    n_child_node_N = 0
    do cc=1,n_child_node_C
        if (adtree(current_node(cc))%nentry .GT. node_minDIVsize) then !Divide 

            !Find node median in target coordinate -> ensures balanced tree
            parent_median = sum(tvtx(adtree(current_node(cc))%entry(:),cindex))/real(adtree(current_node(cc))%nentry,dp)

            !Set median of parent node
            adtree(current_node(cc))%median = parent_median

            !Set coordinate index split of parent
            adtree(current_node(cc))%cdntindex = cindex

            !Initialise child links of this node 
            adtree(current_node(cc))%children(:) = 0

            !Sort into child cells
            child_sort(:) = 0
            Nchild1 = 0
            Nchild2 = 0
            do ii=1,adtree(current_node(cc))%nentry
                if (tvtx(adtree(current_node(cc))%entry(ii),cindex) .GE. parent_median) then !Child 1 (> parent_median)
                    child_sort(ii) = 1
                    Nchild1 = Nchild1 + 1
                else !Child 2 (< parent_median)
                    child_sort(ii) = 2
                    Nchild2 = Nchild2 + 1
                end if
            end do

            !Create reqired child nodes if more than 1 item in current node and they are split
            if (adtree(current_node(cc))%nentry .GT. 1) then
                if((Nchild1 .GT. 1) .AND. (Nchild2 .GT. 1)) then

                    !Child 1 (greater than median) --------------------------
                    allocate(adtree(childins)%entry(Nchild1))
                    adtree(current_node(cc))%children(1) = childins !Set child of parent
                    adtree(childins)%nentry = Nchild1 !Set number of entries of child
                
                    !Insert entries
                    valins = 1
                    do ii=1,adtree(current_node(cc))%nentry
                        if (child_sort(ii) == 1) then
                            adtree(childins)%entry(valins) = adtree(current_node(cc))%entry(ii)
                            valins = valins + 1
                        end if
                    end do 

                    !Set dimension bounds
                    allocate(adtree(childins)%ai(ndim))
                    allocate(adtree(childins)%bi(ndim))
                    do ii=1,ndim
                        adtree(childins)%ai(ii) = minval(tvtx(adtree(childins)%entry(:),ii))
                        adtree(childins)%bi(ii) = maxval(tvtx(adtree(childins)%entry(:),ii))
                    end do
                    
                    !Increment child count
                    new_node(n_child_node_N+1) = childins
                    childins = childins + 1
                    n_child_node_N = n_child_node_N + 1


                    !Child 2 (less than median) --------------------------
                    allocate(adtree(childins)%entry(Nchild2))
                    adtree(current_node(cc))%children(2) = childins !Set child of parent
                    adtree(childins)%nentry = Nchild2 !Set number of entries of child

                    !Insert entries
                    valins = 1
                    do ii=1,adtree(current_node(cc))%nentry
                        if (child_sort(ii) == 2) then
                            adtree(childins)%entry(valins) = adtree(current_node(cc))%entry(ii)
                            valins = valins + 1
                        end if
                    end do 
                    
                    !Set dimension bounds
                    allocate(adtree(childins)%ai(ndim))
                    allocate(adtree(childins)%bi(ndim))
                    do ii=1,ndim
                        adtree(childins)%ai(ii) = minval(tvtx(adtree(childins)%entry(:),ii))
                        adtree(childins)%bi(ii) = maxval(tvtx(adtree(childins)%entry(:),ii))
                    end do

                    !Increment child count
                    new_node(n_child_node_N+1) = childins
                    childins = childins + 1
                    n_child_node_N = n_child_node_N + 1
                end if
            end if
        else !Set child links to zero 
            adtree(current_node(cc))%children(:) = 0
        end if 
    end do

    !Exit if no new children are constructed 
    if (n_child_node_N == 0) then
        exit 
    end if

    !Display
    if (disp_toggle == 1) then
        write(*,'(A,I3,A,I2,A)') '    level :', rr ,' {',cindex,'}'
    end if 
    
    !Update current child search index list
    current_node(1:n_child_node_N) = new_node(1:n_child_node_N)
    n_child_node_C = n_child_node_N

    !Increment depth
    depth = depth + 1
end do
Nnode = childins-1

!Construct final tree 
allocate(adtreeI(Nnode))
do rr=1,Nnode

    !Transfer data 
    allocate(adtreeI(rr)%entry(adtree(rr)%nentry))
    allocate(adtreeI(rr)%ai(ndim))
    allocate(adtreeI(rr)%bi(ndim))
    adtreeI(rr)%nentry = adtree(rr)%nentry
    adtreeI(rr)%entry(:) = adtree(rr)%entry(:)
    adtreeI(rr)%ai(:) = adtree(rr)%ai(:)
    adtreeI(rr)%bi(:) = adtree(rr)%bi(:)
    adtreeI(rr)%children(:) = adtree(rr)%children(:)
    adtreeI(rr)%median = adtree(rr)%median 
    adtreeI(rr)%cdntindex = adtree(rr)%cdntindex
end do 

!Exit display
if (disp_toggle == 1) then
    write(*,'(A,I6,A)') '    {constructed ',childins-1,' nodes}'
end if 
return
end subroutine ad_tree_construct




!Subroutine To Search adtree For Intersection Candidates (4D) =================================================
subroutine ad_tree_search4(node_select,nselected,adtree,ndim,depth_search,tgt_ab,Nnode)
implicit none

!Variables - Import
integer(in) :: depth_search,ndim,nselected,Nnode
integer(in), dimension(:) :: node_select
real(dp), dimension(:,:) :: tgt_ab
type(ADtree_data), dimension(:) :: adtree 

!Variables - Local 
integer(in) :: rr,cc,n_child_node_C,n_child_node_N,node,Nvisit
integer(in) :: current_child(Nnode),new_child(Nnode)
integer(in) :: node_overlap

!Initialise
nselected = 0
Nvisit = 0
!node_select(:) = 0

!Set current node to search
current_child(1) = 1
n_child_node_C = 1

!Search tree
do rr=1,depth_search*ndim

    !Reset new number of chuld nodes to search
    n_child_node_N = 0

    !Check each targeted node
    do cc=1,n_child_node_C

        !Index of node to search
        node = current_child(cc)
        Nvisit = Nvisit + 1

        !Check if node zone overlaps with target zone 
        node_overlap = 0 !Default no overlap 
        if ((adtree(node)%ai(1) .LE. tgt_ab(1,2)) .AND. (adtree(node)%ai(2) .LE. tgt_ab(2,2))) then !node min LEQ bb max
            if ((adtree(node)%bi(3) .GE. tgt_ab(3,1)) .AND. (adtree(node)%bi(4) .GE. tgt_ab(4,1))) then !node max GEQ bb min
                node_overlap = 1
            end if
        end if

        !Cases
        if (node_overlap == 1) then !Node and target zone overlap -> if node has no children then select this node to output else search its children
            if ((adtree(node)%children(1) == 0) .AND. (adtree(node)%children(2) == 0)) then !No children so output node
                node_select(nselected+1) = node
                nselected = nselected + 1 
            else !Check if node has children to search 
                if (adtree(node)%children(1) .NE. 0) then 
                    new_child(n_child_node_N + 1) = adtree(node)%children(1)
                    n_child_node_N = n_child_node_N + 1
                end if
                if (adtree(node)%children(2) .NE. 0) then 
                    new_child(n_child_node_N + 1) = adtree(node)%children(2)
                    n_child_node_N = n_child_node_N + 1
                end if
            end if  
        end if
    end do

    !Exit if no new children are selected to be searched 
    if (n_child_node_N == 0) then
        exit 
    end if

    !Update current child node check index list
    current_child(1:n_child_node_N) = new_child(1:n_child_node_N)
    n_child_node_C = n_child_node_N
end do 
return
end subroutine ad_tree_search4




!Subroutine To Search adtree For Intersection Candidates (6D) =================================================
subroutine ad_tree_search6(node_select,nselected,adtree,ndim,depth_search,tgt_ab,Nnode)
implicit none

!Variables - Import
integer(in) :: depth_search,ndim,nselected,Nnode
integer(in), dimension(:) :: node_select
real(dp), dimension(:,:) :: tgt_ab
type(ADtree_data), dimension(:) :: adtree 

!Variables - Local 
integer(in) :: rr,cc,n_child_node_C,n_child_node_N,node,Nvisit
integer(in) :: current_child(Nnode),new_child(Nnode)
integer(in) :: node_overlap

!Initialise
nselected = 0
Nvisit = 0
!node_select(:) = 0

!Set current node to search
current_child(1) = 1
n_child_node_C = 1

!Search tree
do rr=1,depth_search*ndim

    !Reset new number of chuld nodes to search
    n_child_node_N = 0

    !Check each targeted node
    do cc=1,n_child_node_C

        !Index of node to search
        node = current_child(cc)
        Nvisit = Nvisit + 1

        !Check if node zone overlaps with target zone 
        node_overlap = 0 !Default no overlap 
        if ((adtree(node)%ai(1) .LE. tgt_ab(1,2)) .AND. (adtree(node)%ai(2) .LE. tgt_ab(2,2)) &
        .AND. (adtree(node)%ai(3) .LE. tgt_ab(3,2))) then !node min LEQ bb max
            if ((adtree(node)%bi(4) .GE. tgt_ab(4,1)) .AND. (adtree(node)%bi(5) .GE. tgt_ab(5,1)) &
            .AND. (adtree(node)%bi(6) .GE. tgt_ab(6,1))) then !node max GEQ bb min
                node_overlap = 1
            end if
        end if

        !Cases
        if (node_overlap == 1) then !Node and target zone overlap -> if node has no children then select this node to output else search its children
            if ((adtree(node)%children(1) == 0) .AND. (adtree(node)%children(2) == 0)) then !No children so output node
                node_select(nselected+1) = node
                nselected = nselected + 1 
            else !Check if node has children to search 
                if (adtree(node)%children(1) .NE. 0) then 
                    new_child(n_child_node_N + 1) = adtree(node)%children(1)
                    n_child_node_N = n_child_node_N + 1
                end if
                if (adtree(node)%children(2) .NE. 0) then 
                    new_child(n_child_node_N + 1) = adtree(node)%children(2)
                    n_child_node_N = n_child_node_N + 1
                end if
            end if  
        end if
    end do

    !Exit if no new children are selected to be searched 
    if (n_child_node_N == 0) then
        exit 
    end if

    !Update current child node check index list
    current_child(1:n_child_node_N) = new_child(1:n_child_node_N)
    n_child_node_C = n_child_node_N
end do 
return
end subroutine ad_tree_search6




!Subroutine To Search adtree For Existance Of Any Intersection Candidates =================================================
subroutine ad_tree_search_existance4(exist,adtree,ndim,depth_search,tgt_ab,Nnode)
implicit none

!Variables - Import
integer(in) :: depth_search,ndim,exist,Nnode
real(dp), dimension(:,:) :: tgt_ab
type(ADtree_data), dimension(:) :: adtree 

!Variables - Local 
integer(in) :: rr,cc,n_child_node_C,n_child_node_N,node,Nvisit
integer(in) :: current_child(Nnode),new_child(Nnode)
integer(in) :: node_overlap

!Initialise
exist = 0
Nvisit = 0

!Set current child to search
current_child(1) = 1
n_child_node_C = 1

!Search tree
do rr=1,depth_search*ndim

    !Reset new number of chuld nodes to search
    n_child_node_N = 0

    !Check each targeted node
    do cc=1,n_child_node_C

        !Index of node to search
        node = current_child(cc)
        Nvisit = Nvisit + 1

        !Check if node zone overlaps with target zone 
        node_overlap = 0 !Default no overlap 
        if ((adtree(node)%ai(1) .LE. tgt_ab(1,2)) .AND. (adtree(node)%ai(2) .LE. tgt_ab(2,2))) then !node min LEQ bb max
            if ((adtree(node)%bi(3) .GE. tgt_ab(3,1)) .AND. (adtree(node)%bi(4) .GE. tgt_ab(4,1))) then !node max GEQ bb min
                node_overlap = 1
            end if
        end if

        !Cases
        if (node_overlap == 1) then !Node and target zone overlap -> if node has no children then select this node to output else search its children
            if ((adtree(node)%children(1) == 0) .AND. (adtree(node)%children(2) == 0)) then !No children so possible candidate so return 
                exist = 1
                exit
            else !Check if node has children to search 
                if (adtree(node)%children(1) .NE. 0) then 
                    new_child(n_child_node_N + 1) = adtree(node)%children(1)
                    n_child_node_N = n_child_node_N + 1
                end if
                if (adtree(node)%children(2) .NE. 0) then 
                    new_child(n_child_node_N + 1) = adtree(node)%children(2)
                    n_child_node_N = n_child_node_N + 1
                end if
            end if  
        end if
    end do

    !Exit if no new children are selected to be searched 
    if (n_child_node_N == 0) then
        exit 
    end if

    !Update current child node check index list
    current_child(1:n_child_node_N) = new_child(1:n_child_node_N)
    n_child_node_C = n_child_node_N
end do
return
end subroutine ad_tree_search_existance4




!Subroutine to search adtree for node that is closest to target point =================================================
subroutine ad_tree_search_closest_nodes_2point4(Vtest,distref,node_select,nselected,adtree,ndim,depth_search,tgt_ab,Nnode)
implicit none

!Variables - Import
integer(in) :: depth_search,ndim,nselected,Nnode
integer(in), dimension(:) :: node_select
real(dp) :: distref
real(dp) :: Vtest(2)
real(dp), dimension(:,:) :: tgt_ab
type(ADtree_data), dimension(:) :: adtree 

!Variables - Local 
integer(in) :: rr,cc,cn,n_child_node_C,n_child_node_N,node,Nvisit
integer(in) :: current_child(Nnode),new_child(Nnode)
integer(in) :: node_overlap,nodeminD
real(dp) :: distmin,distchild 
real(dp) :: v1(2),v2(2),v3(2),v4(2),child_dists(2)

!Initialise
nselected = 0
Nvisit = 0
node_select(:) = 0

!Set initial distance 
distmin = distref

!Set current node to search
current_child(1) = 1
n_child_node_C = 1
nodeminD = 1

!Search tree
do rr=1,depth_search*ndim

    !Reset new number of chuld nodes to search
    n_child_node_N = 0

    !Check each targeted node
    do cc=1,n_child_node_C

        !Index of node to search
        node = current_child(cc)
        Nvisit = Nvisit + 1

        !Check if node zone overlaps with target zone 
        node_overlap = 0 !Default no overlap 
        if ((adtree(node)%ai(1) .LE. tgt_ab(1,2)) .AND. (adtree(node)%ai(2) .LE. tgt_ab(2,2))) then !node min LEQ bb max
            if ((adtree(node)%bi(3) .GE. tgt_ab(3,1)) .AND. (adtree(node)%bi(4) .GE. tgt_ab(4,1))) then !node max GEQ bb min
                node_overlap = 1
            end if
        end if

        !Cases
        if (node_overlap == 1) then !Node and target zone overlap -> if node has no children then select this node to output else search its children
            if ((adtree(node)%children(1) == 0) .AND. (adtree(node)%children(2) == 0)) then !No children so add node to output list 
                node_select(nselected+1) = node
                nselected = nselected + 1 
            else !Check if node has children to search 
                if (adtree(node)%children(1) .NE. 0) then 
                    new_child(n_child_node_N + 1) = adtree(node)%children(1)
                    n_child_node_N = n_child_node_N + 1
                end if
                if (adtree(node)%children(2) .NE. 0) then 
                    new_child(n_child_node_N + 1) = adtree(node)%children(2)
                    n_child_node_N = n_child_node_N + 1
                end if
            end if  
        else !Node and target zone do not overlap -> if node has children select closest child to search if closer than minimum current distance to point distmin
            if ((adtree(node)%children(1) == 0) .AND. (adtree(node)%children(2) == 0)) then !No children -> store this node as the minimum distance node so far if closest to point 

                !Node bounding box 
                v1(1) = adtree(node)%ai(1) !xmin 
                v1(2) = adtree(node)%ai(2) !ymin 
                v2(1) = adtree(node)%ai(1) !xmin
                v2(2) = adtree(node)%bi(4) !ymax
                v3(1) = adtree(node)%bi(3) !xmax
                v3(2) = adtree(node)%bi(4) !ymax
                v4(1) = adtree(node)%bi(3) !xmax
                v4(2) = adtree(node)%ai(2) !ymin

                !Distance to point 
                distchild = mindist_point2rectad(Vtest,v1,v2,v3,v4)

                !Test to select
                if (distchild .LT. distmin) then 
                    nodeminD = node
                    distmin = distchild
                end if
            else !Check if node has children to search -> select child closest to point to search 
                do cn=1,2
                    if (adtree(node)%children(cn) .NE. 0) then 

                        !Node bounding box 
                        v1(1) = adtree(adtree(node)%children(cn))%ai(1) !xmin 
                        v1(2) = adtree(adtree(node)%children(cn))%ai(2) !ymin 
                        v2(1) = adtree(adtree(node)%children(cn))%ai(1) !xmin
                        v2(2) = adtree(adtree(node)%children(cn))%bi(4) !ymax
                        v3(1) = adtree(adtree(node)%children(cn))%bi(3) !xmax
                        v3(2) = adtree(adtree(node)%children(cn))%bi(4) !ymax
                        v4(1) = adtree(adtree(node)%children(cn))%bi(3) !xmax
                        v4(2) = adtree(adtree(node)%children(cn))%ai(2) !ymin

                        !Distance to point 
                        child_dists(cn) = mindist_point2rectad(Vtest,v1,v2,v3,v4)
                    end if
                end do 
                if (child_dists(1) .LT. child_dists(2)) then 
                    new_child(n_child_node_N + 1) = adtree(node)%children(1)
                    n_child_node_N = n_child_node_N + 1
                elseif (child_dists(2) .LT. child_dists(1)) then 
                    new_child(n_child_node_N + 1) = adtree(node)%children(2)
                    n_child_node_N = n_child_node_N + 1
                else
                    new_child(n_child_node_N + 1) = adtree(node)%children(1)
                    n_child_node_N = n_child_node_N + 1
                    new_child(n_child_node_N + 1) = adtree(node)%children(2)
                    n_child_node_N = n_child_node_N + 1
                end if 

            end if 
        end if
    end do 

    !Exit if no new children are selected to be searched 
    if (n_child_node_N == 0) then
        exit 
    end if

    !Update current child node check index list
    current_child(1:n_child_node_N) = new_child(1:n_child_node_N)
    n_child_node_C = n_child_node_N
end do 

!Add minimum distance node to node_select
node_select(nselected+1) = nodeminD
nselected = nselected + 1
return 
end subroutine ad_tree_search_closest_nodes_2point4




!Subroutine to search adtree for node that is closest to target zone =================================================
subroutine ad_tree_search_closest_nodes_2zone4(vz1,vz2,vz3,vz4,distref,node_select,nselected,adtree,ndim,depth_search,Nnode)
implicit none

!Variables - Import
integer(in) :: depth_search,ndim,nselected,Nnode
integer(in), dimension(:) :: node_select
real(dp) :: distref
real(dp) :: vz1(2),vz2(2),vz3(2),vz4(2)
type(ADtree_data), dimension(:) :: adtree 

!Variables - Local 
integer(in) :: rr,cc,cn,n_child_node_C,n_child_node_N,node,Nvisit
integer(in) :: current_child(Nnode),new_child(Nnode)
integer(in) :: node_overlap,nodeminD
real(dp) :: distmin,distnode,Zxmin,Zxmax,Zymin,Zymax
real(dp) :: v1(2),v2(2),v3(2),v4(2),child_dists(2)

!Initialise
nselected = 0
Nvisit = 0
node_select(:) = 0

!Set initial distance 
distmin = distref

!Set current node to search
current_child(1) = 1
n_child_node_C = 1
nodeminD = 1

!Set zone bounds
Zxmin = min(vz1(1),vz2(1),vz3(1),vz4(1))
Zxmax = max(vz1(1),vz2(1),vz3(1),vz4(1))
Zymin = min(vz1(2),vz2(2),vz3(2),vz4(2))
Zymax = max(vz1(2),vz2(2),vz3(2),vz4(2))

!Search tree
do rr=1,depth_search*ndim

    !Reset new number of chuld nodes to search
    n_child_node_N = 0

    !Check each targeted node
    do cc=1,n_child_node_C

        !Index of node to search
        node = current_child(cc)
        Nvisit = Nvisit + 1

        !Node bounding box 
        v1(1) = adtree(node)%ai(1) !xmin 
        v1(2) = adtree(node)%ai(2) !ymin 
        v2(1) = adtree(node)%ai(1) !xmin
        v2(2) = adtree(node)%bi(4) !ymax
        v3(1) = adtree(node)%bi(3) !xmax
        v3(2) = adtree(node)%bi(4) !ymax
        v4(1) = adtree(node)%bi(3) !xmax
        v4(2) = adtree(node)%ai(2) !ymin

        !Check if node zone overlaps with target zone 
        node_overlap = 0 !Default no overlap 
        if ((v1(1) .LE. Zxmax) .AND. (v1(2) .LE. Zymax)) then !node min LEQ zone max
            if ((v3(1) .GE. Zxmin) .AND. (v3(2) .GE. Zymin)) then !node max GEQ zone min
                node_overlap = 1
            end if
        end if

        !Cases
        if (node_overlap == 1) then !Node and target zone overlap -> if node has no children then select this node to output else search its children
            if ((adtree(node)%children(1) == 0) .AND. (adtree(node)%children(2) == 0)) then !No children so add node to output list 
                node_select(nselected+1) = node
                nselected = nselected + 1 
            else !Check if node has children to search 
                if (adtree(node)%children(1) .NE. 0) then 
                    new_child(n_child_node_N + 1) = adtree(node)%children(1)
                    n_child_node_N = n_child_node_N + 1
                end if
                if (adtree(node)%children(2) .NE. 0) then 
                    new_child(n_child_node_N + 1) = adtree(node)%children(2)
                    n_child_node_N = n_child_node_N + 1
                end if
            end if  
        else !Node and target zone do not overlap -> select new nodes to search 
            if ((adtree(node)%children(1) == 0) .AND. (adtree(node)%children(2) == 0)) then !No children -> store this node as the minimum distance node so far if closest to point 

                !Distance between node and zone 
                distnode = mindist_rec2recad(v1,v2,v3,v4,vz1,vz2,vz3,vz4)
   
                !Test to select
                if (distnode .LT. distmin) then 
                    nodeminD = node
                    distmin = distnode
                end if
            else !Check if node has children to search -> select child closest to point to search 
                do cn=1,2
                    if (adtree(node)%children(cn) .NE. 0) then 

                        !Node bounding box 
                        v1(1) = adtree(adtree(node)%children(cn))%ai(1) !xmin 
                        v1(2) = adtree(adtree(node)%children(cn))%ai(2) !ymin 
                        v2(1) = adtree(adtree(node)%children(cn))%ai(1) !xmin
                        v2(2) = adtree(adtree(node)%children(cn))%bi(4) !ymax
                        v3(1) = adtree(adtree(node)%children(cn))%bi(3) !xmax
                        v3(2) = adtree(adtree(node)%children(cn))%bi(4) !ymax
                        v4(1) = adtree(adtree(node)%children(cn))%bi(3) !xmax
                        v4(2) = adtree(adtree(node)%children(cn))%ai(2) !ymin

                        !Distance between node and zone 
                        child_dists(cn) = mindist_rec2recad(v1,v2,v3,v4,vz1,vz2,vz3,vz4)
                    end if
                end do 
                if (child_dists(1) .LT. child_dists(2)) then 
                    new_child(n_child_node_N + 1) = adtree(node)%children(1)
                    n_child_node_N = n_child_node_N + 1
                elseif (child_dists(2) .LT. child_dists(1)) then 
                    new_child(n_child_node_N + 1) = adtree(node)%children(2)
                    n_child_node_N = n_child_node_N + 1
                else
                    new_child(n_child_node_N + 1) = adtree(node)%children(1)
                    n_child_node_N = n_child_node_N + 1
                    new_child(n_child_node_N + 1) = adtree(node)%children(2)
                    n_child_node_N = n_child_node_N + 1
                end if 

            end if 
        end if
    end do 

    !Exit if no new children are selected to be searched 
    if (n_child_node_N == 0) then
        exit 
    end if

    !Update current child node check index list
    current_child(1:n_child_node_N) = new_child(1:n_child_node_N)
    n_child_node_C = n_child_node_N
end do 

!Add minimum distance node to node_select
node_select(nselected+1) = nodeminD
nselected = nselected + 1
return 
end subroutine ad_tree_search_closest_nodes_2zone4




!Subroutine To Search adtree For Intersection Candidates With Arbitrary Rectangle =================================================
subroutine ad_tree_search_rectangle_ovlp4(node_select,nselected,rv1,rv2,rv3,rv4,adtree,ndim,depth_search,Nnode)
implicit none

!Variables - Import
integer(in) :: depth_search,ndim,nselected,Nnode
integer(in), dimension(:) :: node_select
real(dp) :: rv1(2),rv2(2),rv3(2),rv4(2)
type(ADtree_data), dimension(:) :: adtree 

!Variables - Local 
integer(in) :: rr,cc,n_child_node_C,n_child_node_N,node,Nvisit
integer(in) :: current_child(Nnode),new_child(Nnode)
integer(in) :: node_overlap,olstatus
real(dp) :: v1(2),v2(2),v3(2),v4(2)

!Initialise
nselected = 0
Nvisit = 0
node_select(:) = 0

!Set current node to search
current_child(1) = 1
n_child_node_C = 1

!Search tree
do rr=1,depth_search*ndim

    !Reset new number of chuld nodes to search
    n_child_node_N = 0

    !Check each targeted node
    do cc=1,n_child_node_C

        !Index of node to search
        node = current_child(cc)
        Nvisit = Nvisit + 1

        !Node bounding box
        v1(1) = adtree(node)%ai(1) !xmin 
        v1(2) = adtree(node)%ai(2) !ymin 
        v2(1) = adtree(node)%ai(1) !xmin
        v2(2) = adtree(node)%bi(4) !ymax
        v3(1) = adtree(node)%bi(3) !xmax
        v3(2) = adtree(node)%bi(4) !ymax
        v4(1) = adtree(node)%bi(3) !xmax
        v4(2) = adtree(node)%ai(2) !ymin

        !Check if node zone overlaps with target rectangle 
        node_overlap = 0 !Default no overlap 
        olstatus = rec_zone_ovlp_bool_2d(v1,v2,v3,v4,rv1,rv2,rv3,rv4)
        if (olstatus .NE. 0) then 
            node_overlap = 1
        end if 

        !Cases
        if (node_overlap == 1) then !Node and target zone overlap -> if node has no children then select this node to output else search its children
            if ((adtree(node)%children(1) == 0) .AND. (adtree(node)%children(2) == 0)) then !No children so output node
                node_select(nselected+1) = node
                nselected = nselected + 1 
            else !Check if node has children to search 
                if (adtree(node)%children(1) .NE. 0) then 
                    new_child(n_child_node_N + 1) = adtree(node)%children(1)
                    n_child_node_N = n_child_node_N + 1
                end if
                if (adtree(node)%children(2) .NE. 0) then 
                    new_child(n_child_node_N + 1) = adtree(node)%children(2)
                    n_child_node_N = n_child_node_N + 1
                end if
            end if  
        end if
    end do

    !Exit if no new children are selected to be searched 
    if (n_child_node_N == 0) then
        exit 
    end if

    !Update current child node check index list
    current_child(1:n_child_node_N) = new_child(1:n_child_node_N)
    n_child_node_C = n_child_node_N
end do 
return
end subroutine ad_tree_search_rectangle_ovlp4




!================================================================================================
!Geometry functions  ============================================================================
!================================================================================================

!Rectangle zone overlap check 2D
function rec_zone_ovlp_bool_2d(zv1,zv2,zv3,zv4,rv1,rv2,rv3,rv4) result(ovlp)
implicit none 

!Variables - Import
integer(in) :: ovlp
real(dp), dimension(:) :: rv1,rv2,rv3,rv4
real(dp), dimension(:) :: zv1,zv2,zv3,zv4

!Variables - Local 
integer(in) :: seg_status 
real(dp) :: zxmax,zxmin,zymax,zymin
real(dp) :: crel(2)

!Initialise 
ovlp = 0

!Check if any vertex of the rectangle is within the zone 
zxmax = max(zv1(1),zv2(1),zv3(1),zv4(1))
zxmin = min(zv1(1),zv2(1),zv3(1),zv4(1))
zymax = max(zv1(2),zv2(2),zv3(2),zv4(2))
zymin = min(zv1(2),zv2(2),zv3(2),zv4(2))
if ((rv1(1) .GE. zxmin) .AND. (rv1(1) .LE. zxmax)) then 
    if ((rv1(2) .GE. zymin) .AND. (rv1(2) .LE. zymax)) then 
        ovlp = 1
        return 
    end if 
end if 
if ((rv2(1) .GE. zxmin) .AND. (rv2(1) .LE. zxmax)) then 
    if ((rv2(2) .GE. zymin) .AND. (rv2(2) .LE. zymax)) then 
        ovlp = 1
        return
    end if 
end if 
if ((rv3(1) .GE. zxmin) .AND. (rv3(1) .LE. zxmax)) then 
    if ((rv3(2) .GE. zymin) .AND. (rv3(2) .LE. zymax)) then 
        ovlp = 1
        return
    end if 
end if 
if ((rv4(1) .GE. zxmin) .AND. (rv4(1) .LE. zxmax)) then 
    if ((rv4(2) .GE. zymin) .AND. (rv4(2) .LE. zymax)) then 
        ovlp = 1
        return
    end if 
end if 

!Check if any vertex of the zone is within the rectangle 
crel = relative_coords_vinrec_2d(zv1,rv1,rv2,rv4)
if ((crel(1) .GE. 0.0d0) .AND. (crel(1) .LE. 1.0d0)) then 
    if ((crel(2) .GE. 0.0d0) .AND. (crel(2) .LE. 1.0d0)) then 
        ovlp = 2
        return
    end if 
end if 
crel = relative_coords_vinrec_2d(zv2,rv1,rv2,rv4)
if ((crel(1) .GE. 0.0d0) .AND. (crel(1) .LE. 1.0d0)) then 
    if ((crel(2) .GE. 0.0d0) .AND. (crel(2) .LE. 1.0d0)) then 
        ovlp = 2
        return
    end if 
end if 
crel = relative_coords_vinrec_2d(zv3,rv1,rv2,rv4)
if ((crel(1) .GE. 0.0d0) .AND. (crel(1) .LE. 1.0d0)) then 
    if ((crel(2) .GE. 0.0d0) .AND. (crel(2) .LE. 1.0d0)) then 
        ovlp = 2
        return
    end if 
end if 
crel = relative_coords_vinrec_2d(zv4,rv1,rv2,rv4)
if ((crel(1) .GE. 0.0d0) .AND. (crel(1) .LE. 1.0d0)) then 
    if ((crel(2) .GE. 0.0d0) .AND. (crel(2) .LE. 1.0d0)) then 
        ovlp = 2
        return
    end if 
end if 

!Check if any edge of the rectangle intersects any edge of the zone 
!R edge 1
seg_status = seg_seg_intersect_boolad(rv1,rv2,zv1,zv2)
if (seg_status .GT. 0) then 
    ovlp = 3 
    return
end if 
seg_status = seg_seg_intersect_boolad(rv1,rv2,zv2,zv3)
if (seg_status .GT. 0) then 
    ovlp = 3 
    return
end if 
seg_status = seg_seg_intersect_boolad(rv1,rv2,zv3,zv4)
if (seg_status .GT. 0) then 
    ovlp = 3 
    return
end if 
seg_status = seg_seg_intersect_boolad(rv1,rv2,zv4,zv1)
if (seg_status .GT. 0) then 
    ovlp = 3 
    return
end if 

!R edge 2
seg_status = seg_seg_intersect_boolad(rv2,rv3,zv1,zv2)
if (seg_status .GT. 0) then 
    ovlp = 3 
    return
end if 
seg_status = seg_seg_intersect_boolad(rv2,rv3,zv2,zv3)
if (seg_status .GT. 0) then 
    ovlp = 3 
    return
end if 
seg_status = seg_seg_intersect_boolad(rv2,rv3,zv3,zv4)
if (seg_status .GT. 0) then 
    ovlp = 3 
    return
end if 
seg_status = seg_seg_intersect_boolad(rv2,rv3,zv4,zv1)
if (seg_status .GT. 0) then 
    ovlp = 3 
    return
end if 

!R edge 3
seg_status = seg_seg_intersect_boolad(rv3,rv4,zv1,zv2)
if (seg_status .GT. 0) then 
    ovlp = 3 
    return
end if 
seg_status = seg_seg_intersect_boolad(rv3,rv4,zv2,zv3)
if (seg_status .GT. 0) then 
    ovlp = 3 
    return
end if 
seg_status = seg_seg_intersect_boolad(rv3,rv4,zv3,zv4)
if (seg_status .GT. 0) then 
    ovlp = 3 
    return
end if 
seg_status = seg_seg_intersect_boolad(rv3,rv4,zv4,zv1)
if (seg_status .GT. 0) then 
    ovlp = 3 
    return
end if 

!R edge 4
seg_status = seg_seg_intersect_boolad(rv4,rv1,zv1,zv2)
if (seg_status .GT. 0) then 
    ovlp = 3 
    return
end if 
seg_status = seg_seg_intersect_boolad(rv4,rv1,zv2,zv3)
if (seg_status .GT. 0) then 
    ovlp = 3 
    return
end if 
seg_status = seg_seg_intersect_boolad(rv4,rv1,zv3,zv4)
if (seg_status .GT. 0) then 
    ovlp = 3 
    return
end if 
seg_status = seg_seg_intersect_boolad(rv4,rv1,zv4,zv1)
if (seg_status .GT. 0) then 
    ovlp = 3 
    return
end if 
return 
end function rec_zone_ovlp_bool_2d




!Find relative coordinates of a point within a rectangle 2D ===========================
function relative_coords_vinrec_2d(vtgt,rv1,rv2,rv4) result(crel)
implicit none 

!Variables - Import
real(dp) :: crel(2),vtgt(2)
real(dp), dimension(:) :: rv1,rv2,rv4

!Variables - Local 
real(dp) :: normx,normy
real(dp) :: Rxv(2),Ryv(2),vec_ref(2)

!Find relative coordinates 
Rxv(:) = rv2(:) - rv1(:)
Ryv(:) = rv4(:) - rv1(:)
normx = sqrt(Rxv(1)**2 + Rxv(2)**2)
normy = sqrt(Ryv(1)**2 + Ryv(2)**2)
vec_ref(:) = vtgt(:) - rv1(:)
crel(1) = dot_product(vec_ref,Rxv)/(normx*normx)
crel(2) = dot_product(vec_ref,Ryv)/(normy*normy)
return 
end function relative_coords_vinrec_2d




!Build 2D rectangle from centreline ===========================
subroutine rec_from_centrline_2d(rv1,rv2,rv3,rv4,clv1,clv2,recwidth)

!Variables - Import
real(dp) :: recwidth
real(dp), dimension(:) :: clv1,clv2
real(dp), dimension(:), allocatable :: rv1,rv2,rv3,rv4

!Variables - Local 
real(dp) :: dx,dy,lnorm 
real(dp) :: rnorm(2)

!Allocate
allocate(rv1(2))
allocate(rv2(2))
allocate(rv3(2))
allocate(rv4(2))

!Build vertices
dx = clv2(1) - clv1(1)
dy = clv2(2) - clv1(2)
Lnorm = sqrt(dx*dx + dy*dy)
rnorm(1) = dy/Lnorm
rnorm(2) = -dx/Lnorm
rv1(:) = clv1(:) - recwidth*rnorm(:)
rv2(:) = clv1(:) + recwidth*rnorm(:)
rv3(:) = clv2(:) + recwidth*rnorm(:)
rv4(:) = clv2(:) - recwidth*rnorm(:)
end subroutine rec_from_centrline_2d




!Line-Line intersection location calculation -> L1 = v1 v2 | L2 = v3 v4 ===========================
function line_line_intersection_loc_inl1ad(v1,v2,v3,v4) result(vi)
use ieee_arithmetic, only: ieee_value,IEEE_QUIET_NAN
implicit none 

!Variables - Import
real(dp) :: v1(2),v2(2),v3(2),v4(2)

!Variables - Local
real(dp) :: Dval,t

!Result
real(dp) :: vi(2)

!Intersection denominator 
Dval = (v1(1) - v2(1))*(v3(2) - v4(2)) - (v1(2) - v2(2))*(v3(1) - v4(1))

!Test paralell
if (Dval .NE. 0.0d0) then !If not parallel -> find intersection location
    t = ((v1(1) - v3(1))*(v3(2) - v4(2)) - (v1(2) - v3(2))*(v3(1) - v4(1)))/Dval !L1 parameter
    if (t .GT. 1.0d0) then 
        t = 1.0d0 
    elseif (t .LT. 0.0d0) then 
        t = 0.0d0 
    end if
    vi(1) = (v1(1) + t*(v2(1) - v1(1)))
    vi(2) = (v1(2) + t*(v2(2) - v1(2)))
else !Take as not intersecting as parallel  
    vi(:) = ieee_value(1.0d0,IEEE_QUIET_NAN)
end if
return 
end function line_line_intersection_loc_inl1ad
    



!Minimum distance from a point to a line segment function ===========================
function mindist_point2linead(vp,v1,v2) result(vid)
implicit none 

!Variables - Import
real(dp) :: vid(3) !xi | yi | dist 
real(dp) :: vp(2),v1(2),v2(2)

!Variables - Local 
real(dp) :: dx,dy,nx,ny
real(dp) :: vn1(2),vn2(2)

!Find line segment direction 
dx = v2(1) - v1(1)
dy = v2(2) - v1(2)

!Check for zero length segment 
if (dx*dx + dy*dy == 0.0d0) then 

    !Return intersection point is vi = v1 = v2
    vid(1:2) = v1(:)
else

    !Segment normal 
    nx = -dy
    ny = dx

    !Build normal direction vertices 
    vn1(1) = vp(1) + nx 
    vn1(2) = vp(2) + ny 
    vn2(1) = vp(1) - nx 
    vn2(2) = vp(2) - ny

    !Find intersection on base line segment 
    vid(1:2) = line_line_intersection_loc_inl1ad(v1,v2,vn1,vn2) 
end if

!Find distance from this intersect to vp 
vid(3) = sqrt((vp(1) - vid(1))**2 + (vp(2) - vid(2))**2)
return 
end function mindist_point2linead




!Minimum distance from a point to a rectangle function ===========================
function mindist_point2rectad(vp,v1,v2,v3,v4) result(dist)
implicit none 

!Variables - Import
real(dp) :: dist
real(dp) :: vp(2),v1(2),v2(2),v3(2),v4(2)

!Variables - Local 
real(dp) :: diste(4),vide(3)

!Evaluate distance for each edge 
vide = mindist_point2linead(vp,v1,v2)
diste(1) = vide(3)
vide = mindist_point2linead(vp,v2,v3)
diste(2) = vide(3)
vide = mindist_point2linead(vp,v3,v4)
diste(3) = vide(3)
vide = mindist_point2linead(vp,v4,v1)
diste(4) = vide(3)

!select minimum distance 
dist = minval(diste)
return 
end function mindist_point2rectad




!Minimum distance from a line to a line function ===========================
function mindist_line2linead(v1,v2,v3,v4) result(dist)
implicit none 

!Variables - Import
real(dp) :: dist
real(dp) :: v1(2),v2(2),v3(2),v4(2)

!Variables - Local 
integer(in) :: int_type
real(dp) :: dist_e2l(4),vide(3)

!Check if lines intersect or touch 
int_type = seg_seg_intersect_boolad(v1,v2,v3,v4)

!Cases on line-line intersection
if (int_type == 0) then !No intersection -> distance is mimimum from each end of each to the other line
    vide = mindist_point2linead(v1,v3,v4)
    dist_e2l(1) = vide(3) !v1 to v3-v4 line
    vide = mindist_point2linead(v2,v3,v4)
    dist_e2l(2) = vide(3) !v2 to v3-v4 line
    vide = mindist_point2linead(v3,v1,v2)
    dist_e2l(3) = vide(3) !v3 to v1-v2 line
    vide = mindist_point2linead(v4,v1,v2)
    dist_e2l(4) = vide(3) !v4 to v1-v2 line
    dist = minval(dist_e2l)
else !Intersection -> distance is zero 
    dist = 0.0d0 
end if 
return 
end function mindist_line2linead




!Minimum distance from a line to a rectangle function ===========================
function mindist_line2recad(v1,v2,r_1_v1,r_1_v2,r_1_v3,r_1_v4) result(dist)
implicit none 

!Variables - Import
real(dp) :: dist
real(dp) :: r_1_v1(2),r_1_v2(2),r_1_v3(2),r_1_v4(2),v1(2),v2(2)

!Variables - Local 
real(dp) :: dist_edg(4)

!Find distance to each edge from the base line 
dist_edg(1) = mindist_line2linead(v1,v2,r_1_v1,r_1_v2)
dist_edg(2) = mindist_line2linead(v1,v2,r_1_v2,r_1_v3)
dist_edg(3) = mindist_line2linead(v1,v2,r_1_v3,r_1_v4)
dist_edg(4) = mindist_line2linead(v1,v2,r_1_v4,r_1_v1)

!Return minimum distance 
dist = minval(dist_edg)
return 
end function mindist_line2recad




!Minimum distance from a rectangle to a rectangle function ===========================
function mindist_rec2recad(r_1_v1,r_1_v2,r_1_v3,r_1_v4,r_2_v1,r_2_v2,r_2_v3,r_2_v4) result(dist)
implicit none 

!Variables - Import
real(dp) :: dist
real(dp) :: r_1_v1(2),r_1_v2(2),r_1_v3(2),r_1_v4(2)
real(dp) :: r_2_v1(2),r_2_v2(2),r_2_v3(2),r_2_v4(2)

!Variables - Local 
real(dp) :: dist_er1_2_r2(4)

!Find minimum distance from each edge of r1 to rectangle r2   
dist_er1_2_r2(1) = mindist_line2recad(r_1_v1,r_1_v2,r_2_v1,r_2_v2,r_2_v3,r_2_v4)
dist_er1_2_r2(2) = mindist_line2recad(r_1_v2,r_1_v3,r_2_v1,r_2_v2,r_2_v3,r_2_v4)
dist_er1_2_r2(3) = mindist_line2recad(r_1_v3,r_1_v4,r_2_v1,r_2_v2,r_2_v3,r_2_v4)
dist_er1_2_r2(4) = mindist_line2recad(r_1_v4,r_1_v1,r_2_v1,r_2_v2,r_2_v3,r_2_v4)

!Return minimum distance 
dist = minval(dist_er1_2_r2)
return 
end function mindist_rec2recad




!Segment - Segment Intersection Checking Function ===========================
function seg_seg_intersect_boolad(A,B,C,D) result(int_type) !Segments A->B || C->D
implicit none 

!Variables - Import
integer(in) :: int_type
real(dp) :: A(2),B(2),C(2),D(2)

!Variables - Local 
integer(in) :: intA,intB,intC,intD
real(dp) :: cp_abc,cp_abd,cp_cda,cp_cdb
real(dp) :: L1,L2,Afcd,Bfcd,Cfab,Dfab

!Find segment cross products
cp_abc = cp2dad(A,B,C)
cp_abd = cp2dad(A,B,D)
cp_cda = cp2dad(C,D,A)
cp_cdb = cp2dad(C,D,B)

!Classify intersection (0 = none | 1 = proper | 2 = vertex-line touch | 3 = verex-vertex touch | 4 = colinear with overlap)
int_type = 0 
if ((cp_abc*cp_abd .LT. 0.0d0) .AND. (cp_cda*cp_cdb .LT. 0.0d0)) then !Proper intersection (both cross within each others length)
    int_type = 1
elseif ((cp_abc == 0.0d0) .AND. (cp_abd == 0.0d0) .AND. (cp_cda == 0.0d0) .AND. (cp_cdb == 0.0d0)) then !Co-linear
    L1 = norm2(A(:)-B(:))
    L2 = norm2(C(:)-D(:))
    if ((L1 == 0.0d0) .OR. (L2 == 0.0d0)) then !one segment is zero length so ignore this check 
        int_type = 0
    else
        Afcd = dot_product(A(:)-D(:),C(:)-D(:))/(L2**2)
        Bfcd = dot_product(B(:)-D(:),C(:)-D(:))/(L2**2)
        Cfab = dot_product(C(:)-B(:),A(:)-B(:))/(L1**2)
        Dfab = dot_product(D(:)-B(:),A(:)-B(:))/(L1**2)
        if ((Afcd .GT. 0.0d0) .AND. (Afcd .LT. 1.0d0)) then !Contained within segment 
            intA = 4
        elseif ((Afcd == 0.0d0) .OR. (Afcd == 1.0d0)) then !On segment vertex
            intA = 3
        else !Outside segment 
            intA = 0
        end if 
        if ((Bfcd .GT. 0.0d0) .AND. (Bfcd .LT. 1.0d0)) then !Contained within segment 
            intB = 4
        elseif ((Bfcd == 0.0d0) .OR. (Bfcd == 1.0d0)) then !On segment vertex
            intB = 3
        else !Outside segment 
            intB = 0
        end if 
        if ((Cfab .GT. 0.0d0) .AND. (Cfab .LT. 1.0d0)) then !Contained within segment 
            intC = 4
        elseif ((Cfab == 0.0d0) .OR. (Cfab == 1.0d0)) then !On segment vertex
            intC = 3
        else !Outside segment 
            intC = 0
        end if 
        if ((Dfab .GT. 0.0d0) .AND. (Dfab .LT. 1.0d0)) then !Contained within segment 
            intD = 4
        elseif ((Dfab == 0.0d0) .OR. (Dfab == 1.0d0)) then !On segment vertex
            intD = 3
        else !Outside segment 
            intD = 0
        end if 
        if ((intA == 4) .OR. (intB == 4) .OR. (intC == 4) .OR. (intD == 4)) then !Segment - Segment overlap 
            int_type = 4
        elseif ((intA == 3) .OR. (intB == 3) .OR. (intC == 3) .OR. (intD == 3)) then !Vertex - Vertex touch 
            int_type = 3
        else !No intersection
            int_type = 0 
        end if 
    end if 
elseif (cp_abc == 0.0d0) then !C lies on AB
    if (cp_cda*cp_cdb .LT. 0.0d0) then !Within AB
        int_type = 2
    elseif (cp_cda*cp_cdb .GT. 0.0d0) then !Outside AB
        int_type = 0 
    elseif (cp_cda*cp_cdb == 0.0d0) then !On vertex A or B
        int_type = 3
    end if 
elseif (cp_abd == 0.0d0) then !D lies on AB
    if (cp_cda*cp_cdb .LT. 0.0d0) then !Within AB
        int_type = 2
    elseif (cp_cda*cp_cdb .GT. 0.0d0) then !Outside AB
        int_type = 0 
    elseif (cp_cda*cp_cdb == 0.0d0) then !On vertex A or B
        int_type = 3
    end if
elseif (cp_cda == 0.0d0) then !A lies on CD
    if (cp_abc*cp_abd .LT. 0.0d0) then !Within CD
        int_type = 2
    elseif (cp_abc*cp_abd .GT. 0.0d0) then !Outside CD  
        int_type = 0 
    elseif (cp_abc*cp_abd == 0.0d0) then !On vertex C or D
        int_type = 3
    end if
elseif (cp_cdb == 0.0d0) then !B lies on CD
    if (cp_abc*cp_abd .LT. 0.0d0) then !Within CD
        int_type = 2
    elseif (cp_abc*cp_abd .GT. 0.0d0) then !Outside CD  
        int_type = 0 
    elseif (cp_abc*cp_abd == 0.0d0) then !On vertex C or D
        int_type = 3
    end if
else !No intersection within either lines length
    int_type = 0 
end if  
return 
end function seg_seg_intersect_boolad
    



!2D Cross product functions ===========================
function cp2dad(A,B,C) result(val) !vectors (A->B) and (A->C)
implicit none 

!Variables - Import
real(dp) :: val
real(dp) :: A(2),B(2),C(2)

!Evaluate
val = (B(1) - A(1))*(C(2) - A(2)) - (B(2) - A(2))*(C(1) - A(1))
return 
end function cp2dad


end module cellmesh2d_adtree_mod