!2D Geometry Routine Module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 8.4
!Updated 20-02-2024

!Geometry subroutines module
module cellmesh2d_geometry_mod
use cellmesh2d_adtree_mod
use cellmesh2d_utilities_mod

contains


!Set comparison precision function ===========================
function set_cmp_prc(reflen,min_precision_val) result(cval)
implicit none 

!Variables - Import 
real(dp) :: reflen,min_precision_val,cval

!Set comparision zero bound to the maximum of the reference length and the round off error minimum precision bound 
cval = max(reflen,min_precision_val)
return 
end function set_cmp_prc




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




!Segment - Rectangle Intersection Checking Function ===========================
function seg_rectangle_intersect_bool(vp1,vp2,vp3,vp4,A,B) result(int_type)
implicit none 

!Variables - Import
integer(in) :: int_type
real(dp) :: vp1(2),vp2(2),vp3(2),vp4(2),A(2),B(2)

!Variables - Local 
integer(in) :: int_e1,int_e2,int_e3,int_e4,in_A,in_B,int_init
real(dp) :: xmin,xmax,ymin,ymax

!Test intersection of line segment A->B against each edge of the rectangle
int_init = seg_seg_intersect_bool(A,B,vp1,vp2)
if (int_init .NE. 0) then 
    int_e1 = 1
else 
    int_e1 = 0
end if
int_init = seg_seg_intersect_bool(A,B,vp2,vp3)
if (int_init .NE. 0) then 
    int_e2 = 1
else 
    int_e2 = 0
end if
int_init = seg_seg_intersect_bool(A,B,vp3,vp4)
if (int_init .NE. 0) then 
    int_e3 = 1
else 
    int_e3 = 0
end if
int_init = seg_seg_intersect_bool(A,B,vp4,vp1)
if (int_init .NE. 0) then 
    int_e4 = 1
else 
    int_e4 = 0
end if

!Test containment of each segment end 
in_A = 0
in_B = 0 
xmin = min(vp1(1),vp2(1),vp3(1),vp4(1))
xmax = max(vp1(1),vp2(1),vp3(1),vp4(1))
ymin = min(vp1(2),vp2(2),vp3(2),vp4(2))
ymax = max(vp1(2),vp2(2),vp3(2),vp4(2))
if ((A(1) .GE. xmin) .AND. (A(1) .LE. xmax)) then 
    if ((A(2) .GE. ymin) .AND. (A(2) .LE. ymax)) then 
        in_A = 1
    end if    
end if
if ((B(1) .GE. xmin) .AND. (B(1) .LE. xmax)) then 
    if ((B(2) .GE. ymin) .AND. (B(2) .LE. ymax)) then 
        in_B = 1
    end if    
end if

!Set containment status 
if ((int_e1 == 1) .OR. (int_e2 == 1) .OR. (int_e3 == 1) .OR. (int_e4 == 1) .OR. (in_A == 1) .OR. (in_B == 1)) then 
    int_type = 1
else
    int_type = 0 
end if
return 
end function seg_rectangle_intersect_bool




!Segment - Segment Intersection Checking Function ===========================
function seg_seg_intersect_bool(A,B,C,D) result(int_type) !Segments A->B || C->D
implicit none 

!Variables - Import
integer(in) :: int_type
real(dp) :: A(2),B(2),C(2),D(2)

!Variables - Local 
integer(in) :: intA,intB,intC,intD
real(dp) :: cp_abc,cp_abd,cp_cda,cp_cdb
real(dp) :: L1,L2,Afcd,Bfcd,Cfab,Dfab

!Find segment cross products
cp_abc = cp2d(A,B,C)
cp_abd = cp2d(A,B,D)
cp_cda = cp2d(C,D,A)
cp_cdb = cp2d(C,D,B)

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
end function seg_seg_intersect_bool




!2D Cross product functions ===========================
function cp2d(A,B,C) result(val) !vectors (A->B) and (A->C)
implicit none 

!Variables - Import
real(dp) :: val
real(dp) :: A(2),B(2),C(2)

!Evaluate
val = (B(1) - A(1))*(C(2) - A(2)) - (B(2) - A(2))*(C(1) - A(1))
return 
end function cp2d

function cp2dv(vec1N,vec2N) result(val) !vectors (vec1N = (A->B)) and (vec2N = (A->C))
implicit none 

!Variables - Import
real(dp) :: val
real(dp) :: vec1N(2),vec2N(2)

!Evaluate
val = vec1N(1)*vec2N(2) - vec1N(2)*vec2N(1)
return 
end function cp2dv




!Outer angle function ===========================
function outer_angle(v1,vp,v2) result(ang)
implicit none 

!Variables - Import
real(dp) :: ang
real(dp) :: v1(2),v2(2),vp(2)

!Variables - Local 
real(dp) :: cpval,cosang,pi
real(dp) :: ve1(2),ve2(2)

!Define pi
pi = 4.0d0*atan(1.0d0)

!Edge vectors 
ve1(:) = vp(:) - v1(:)
ve2(:) = vp(:) - v2(:)
ve1(:) = ve1(:)/norm2(ve1(:))
ve2(:) = ve2(:)/norm2(ve2(:))

!Check convexity 
cpval = cp2dv(ve1,ve2)

!Construct angle 
cosang = dot_product(ve1,ve2)
if (cosang .GT. 1.0d0) then 
    cosang = 1.0d0 
elseif (cosang .LT. -1.0d0) then 
    cosang = -1.0d0 
end if 
ang = acos(cosang)
ang = ang*(180.0d0/pi)
if (cpval .GT. 0.0d0) then 
    ang = 360.0d0 - ang
end if 
ang = abs(ang)
return 
end function outer_angle




!Line segment area contribution function ===========================
function Asegment(v1,v2) result(Aseg) !+ve area for CCW oriented shapes
implicit none 

!Variables - Import
real(dp) :: v1(2),v2(2)

!Result
real(dp) :: Aseg

!Segment area
Aseg = 0.5d0*(v1(1)*v2(2) - v2(1)*v1(2))
return 
end function Asegment




!Segment normal function ===========================
function Nsegment(v1,v2) result(Nseg) 
implicit none 

!Variables - Import
real(dp) :: v1(2),v2(2)

!Result
real(dp) :: Nseg(2)

!Variables - Local 
real(dp) :: dx,dy 

!Evaluate
dx = v2(1) - v1(1)
dy = v2(2) - v1(2)
Nseg(1) = -dy 
Nseg(2) = dx
return 
end function Nsegment




!Minimum distant point to edge function ===========================
function min_dist_point_to_edge(ve1,ve2,vp) result(vid)
implicit none 

!Variables - Import
real(dp) :: vp(2),ve1(2),ve2(2),vid(3)

!Variables - Local 
real(dp) :: t,f
real(dp) :: edir(2)

!Evaluate closest point 
edir(:) = ve2(:) - ve1(:)
edir(:) = edir(:)/norm2(edir(:))
t = dot_product(vp - ve1,edir)
vid(1:2) = ve1(:) + t*edir(:)
f = dot_product(vid(1:2) - ve1(:),edir(:))/norm2(ve2(:) - ve1(:))
! f = norm2(vid(1:2) - ve1(:))/norm2(ve2(:) - ve1(:))
if (f .GT. 1.0d0) then 
    vid(1:2) = ve2(:)
elseif (f .LT. 0.0d0) then 
    vid(1:2) = ve1(:)
end if 
vid(3) = norm2(vid(1:2) - vp(:))
return 
end function min_dist_point_to_edge
    



!Line-Line intersection location calculation (bounded) -> L1 = v1 v2 | L2 = v3 v4 ===========================
function line_line_intersection_loc_inl1(v1,v2,v3,v4) result(vi)
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
end function line_line_intersection_loc_inl1




!Local radius of curvature function ===========================
function surface_rcurv(Ninterp,interp_stencil,vertices,cm2dopt) result(Rcurv) 
! use ieee_arithmetic
implicit none 
!Evaluates the surface radius of curvature at the central vertex of interp_stencil

!Variables - Import
integer(in) :: Ninterp
integer(in), dimension(:) :: interp_stencil
real(dp) :: Rcurv
real(dp), dimension(:,:) :: vertices
type(cm2d_options) :: cm2dopt

!Variables - Local
integer(in) :: ii,jj
integer(in) :: interp_pnt
integer(in) :: point_list(Ninterp)
real(dp) :: Rs,dx_ds,dy_ds,dx_ds_2,dy_ds_2,denom
real(dp) :: sdelta(Ninterp)
real(dp) :: vtx_l(Ninterp,2),Rd(Ninterp,Ninterp),Rdi(Ninterp,Ninterp)
real(dp) :: gamma_x(Ninterp,1),gamma_y(Ninterp,1),s(Ninterp,1)
! real(dp) :: spi,vpi(2)

!Set interpolation point 
interp_pnt = ((Ninterp - 1)/2) + 1

!Build local coordinate
s(:,1) = 0.0d0 
do ii=2,Ninterp
    s(ii,1) = norm2(vertices(interp_stencil(ii),:) - vertices(interp_stencil(ii-1),:)) + s(ii-1,1)
end do 
do ii=1,Ninterp
    point_list(ii) = ii
end do 

!Check for coincident vertices 
do ii=2,Ninterp
    if (abs(s(ii,1) - s(ii-1,1)) == 0.0d0) then 
        print *, '** detected two co-incident surface vertices ',interp_stencil(ii),interp_stencil(ii-1)
        Rcurv = 0.0d0 
        return 
    end if 
end do 

!Build local vertex vector
do ii=1,Ninterp
    vtx_l(ii,:) = vertices(interp_stencil(ii),:)
end do 

!Build local dependance matrix
Rd(:,:) = 0.0d0 
do ii=1,Ninterp
    do jj=1,Ninterp
        Rd(ii,jj) = abs(s(ii,1) - s(jj,1))
    end do 
end do
Rs = cm2dopt%RBF_rsup*s(Ninterp,1) 
do ii=1,Ninterp
    do jj=1,Ninterp
        Rd(ii,jj) = wendlandC2(Rd(ii,jj),Rs)
    end do 
end do

!Invert dependance matrix 
call matinv(Rd,Rdi,Ninterp)

!Build x coefficients
gamma_x(:,1) = matmul(Rdi,vtx_l(:,1))

!Build y coefficients
gamma_y(:,1) = matmul(Rdi,vtx_l(:,2))

!First gradients in x and y
do ii=1,Ninterp
    sdelta(ii) = wendlandc2_gradient1(s(interp_pnt,1),s(ii,1),Rs)
end do 
dx_ds = sum(sdelta(:)*gamma_x(:,1))
dy_ds = sum(sdelta(:)*gamma_y(:,1))

!Second gradients in x and y
do ii=1,Ninterp
    sdelta(ii) = wendlandc2_gradient2(abs(s(interp_pnt,1)-s(ii,1)),Rs)
end do 
dx_ds_2 = sum(sdelta(:)*gamma_x(:,1))
dy_ds_2 = sum(sdelta(:)*gamma_y(:,1))

!Calculate radius of curvature at v0 
denom = abs(dx_ds*dy_ds_2 - dy_ds*dx_ds_2) 
if (denom == 0.0d0) then !strait line -> infinite radius 
    Rcurv = ieee_value(1.0d0,IEEE_POSITIVE_INF)
else !curve
    Rcurv = ((dx_ds**2 + dy_ds**2)**1.5d0)/denom
end if 
! print *, 'Rcan = ',Rcurv

!Test export debug -----
! if (interp_stencil(interp_pnt) == 120) then 
!     print *,'Rcurv = ',Rcurv
!     open(11,file='io/interp_curve.dat')
!     do ii=1,1000 
!         spi = ((s(Ninterp) - s(1))/999)*real(ii-1,dp) + s(1)
!         do jj=1,Ninterp
!             sdelta(jj) = wendlandC2(abs(spi - s(jj)),Rs)
!         end do
!         vpi(1) = sum(sdelta(:)*gamma_x(:,1))
!         vpi(2) = sum(sdelta(:)*gamma_y(:,1))
!         write(11,*) vpi(:)
!     end do 
!     close(11)
!     open(11,file='io/interp_curve_vbase.dat')
!         write(11,*) interp_stencil
!     close(11)
! end if 
return 
end function surface_rcurv




!Check geometry for self intersections function ===========================
function is_self_intersecting(surface_mesh,cm2dopt) result(is_selfintersecting)
implicit none 

!Variables - Import
logical :: is_selfintersecting
type(surface_data) :: surface_mesh
type(cm2d_options) :: cm2dopt

!Variables - Local 
integer(in) :: ii,nn,kk
integer(in) :: Ndim,node_minDIVsize,fp,fc,fn,nselected,segval,etgt
integer(in), dimension(:), allocatable:: node_select
real(dp) :: global_target_pad,cpadSZ,zxmin,zxmax,zymin,zymax,zzmin,zzmax
real(dp) :: v1(2),v2(2),vt1(2),vt2(2)
real(dp), dimension(:,:), allocatable :: tvtx
type(tree_data) :: surface_adtree

!Initialise intersection state
is_selfintersecting = .false.

!Set global object bounding box padding for adtree node containement
global_target_pad = 0.0d0

!Adtree number of dimensions (4)
Ndim = 4

!Set minimum divisible node size within the AD tree 
node_minDIVsize = cm2dopt%ADTminNodedivsize

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

!Check each segment for self intersections with all but the two directly adjacent segments 
allocate(node_select(surface_adtree%nnode))
node_select(:) = 0 
nselected = 0 
do ii=1,surface_mesh%nfcs

    !Faces to exclude (previous / current / next)
    fp = surface_mesh%v2f(surface_mesh%faces(ii,1),1)
    fc = ii 
    fn = surface_mesh%v2f(surface_mesh%faces(ii,2),2)

    !Edge ends
    v1(:) = surface_mesh%vertices(surface_mesh%faces(ii,1),:)
    v2(:) = surface_mesh%vertices(surface_mesh%faces(ii,2),:)

    !Set padding size
    cpadSZ = norm2(v2 - v1)

    !Intersection bounding box
    zxmin = min(v1(1),v2(1)) - cpadSZ*cm2dopt%ADTpadding !tgt bounding box -> xmin
    zxmax = max(v1(1),v2(1)) + cpadSZ*cm2dopt%ADTpadding !tgt bounding box -> xmax
    zymin = min(v1(2),v2(2)) - cpadSZ*cm2dopt%ADTpadding !tgt bounding box -> ymin
    zymax = max(v1(2),v2(2)) + cpadSZ*cm2dopt%ADTpadding !tgt bounding box -> ymax

    !Identify any edge bounding boxes that may overlap the edge 
    call search_ADtree(nselected,node_select,surface_adtree,zxmin,zxmax,zymin,zymax,zzmin,zzmax)

    !Check for intersection with all edges selected 
    do nn=1,nselected
        do kk=1,surface_adtree%tree(node_select(nn))%nentry

            !Selected edge
            etgt = surface_adtree%tree(node_select(nn))%entry(kk)

            !If this edge is not any of the excluded three
            if ((etgt .NE. fp) .AND. (etgt .NE. fc) .AND. (etgt .NE. fn)) then 

                !Verticies of the ends of the surface segment 
                vt1(:) = surface_mesh%vertices(surface_mesh%faces(etgt,1),:)
                vt2(:) = surface_mesh%vertices(surface_mesh%faces(etgt,2),:)

                !Intersection type if any 
                segval = seg_seg_intersect_bool(v1,v2,vt1,vt2)

                !If an intersection set tag and exit 
                if (segval .NE. 0) then 
                    is_selfintersecting = .true.
                    exit 
                end if 
            end if 
        end do 
    end do 

    !Exit if self intersection found 
    if (is_selfintersecting) then 
        exit 
    end if 
end do 
return 
end function is_self_intersecting


end module cellmesh2d_geometry_mod