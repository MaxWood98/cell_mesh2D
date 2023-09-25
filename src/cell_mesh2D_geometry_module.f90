!2D Geometry Routine Module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 8.1
!Updated 18-09-2023

!Geometry subroutines module
module cellmesh2d_geometry_mod
use cellmesh2d_adtree_mod
use cellmesh2d_utilities_mod
use ieee_arithmetic !, only: ieee_value,IEEE_QUIET_NAN
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




!Wendland C2 function ===========================
function wendlandc2(d,Rs) result(W)
implicit none 

!Variables - Import
real(dp) :: d,Rs,W

!Set function value 
d = d/Rs 
if (d .GE. 1.0d0) then 
    W = 0.0d0 
else
    W = ((1.0d0 - d)**4)*(4.0d0*d + 1.0d0)
end if 
return 
end function wendlandc2




!Wendland C2 first gradient function ===========================
function wendlandc2_gradient1(s_p,s_b,Rs) result(W)
implicit none 

!Variables - Import
real(dp) :: s_p,s_b,Rs,W

!Variables - Local
real(dp) :: d

!Set function value 
d = abs(s_p - s_b)/Rs 
if (d .GE. 1.0d0) then 
    W = 0.0d0 
else
    W = (4.0d0*(d - 1.0d0)**4 + 4.0d0*(4.0d0*d + 1.0d0)*((d - 1.0d0)**3))*sign(1.0d0,s_p-s_b)
end if 
return 
end function wendlandc2_gradient1




!Wendland C2 second gradient function ===========================
function wendlandc2_gradient2(d,Rs) result(W)
implicit none 

!Variables - Import
real(dp) :: d,Rs,W

!Set function value 
d = d/Rs 
if (d .GE. 1.0d0) then 
    W = 0.0d0 
else
    W = 32.0d0*((d - 1.0d0)**3) + 12.0d0*(4.0d0*d + 1.0d0)*((d - 1.0d0)**2)
end if 
return 
end function wendlandc2_gradient2




!Local radius of curvature function ===========================
function surface_rcurv(Ninterp,interp_stencil,vertices) result(Rcurv) 
use ieee_arithmetic
implicit none 
!Evaluates the surface radius of curvature at the central vertex of interp_stencil

!Variables - Import
integer(in) :: Ninterp
integer(in), dimension(:) :: interp_stencil
real(dp) :: Rcurv
real(dp), dimension(:,:) :: vertices

!Variables - Local
integer(in) :: ii,jj
integer(in) :: interp_pnt
real(dp) :: Rs,dx_ds,dy_ds,dx_ds_2,dy_ds_2,denom
real(dp) :: s(Ninterp),sdelta(Ninterp)
real(dp) :: vtx_l(Ninterp,2),Rd(Ninterp,Ninterp),Rdi(Ninterp,Ninterp)
real(dp) :: gamma_x(Ninterp,1),gamma_y(Ninterp,1)
! real(dp) :: spi,vpi(2)

!Set interpolation point 
interp_pnt = ((Ninterp - 1)/2) + 1

!Build local coordinate
s(:) = 0.0d0 
do ii=2,Ninterp
    s(ii) = norm2(vertices(interp_stencil(ii),:) - vertices(interp_stencil(ii-1),:)) + s(ii-1)
end do 

!Check for coincident vertices 
do ii=2,Ninterp
    if (abs(s(ii) - s(ii-1)) == 0.0d0) then 
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
        Rd(ii,jj) = abs(s(ii) - s(jj))
    end do 
end do
Rs = 80.0d0*s(Ninterp) 
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
    sdelta(ii) = wendlandc2_gradient1(s(interp_pnt),s(ii),Rs)
end do 
dx_ds = sum(sdelta(:)*gamma_x(:,1))
dy_ds = sum(sdelta(:)*gamma_y(:,1))

!Second gradients in x and y
do ii=1,Ninterp
    sdelta(ii) = wendlandc2_gradient2(abs(s(interp_pnt)-s(ii)),Rs)
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


end module cellmesh2d_geometry_mod