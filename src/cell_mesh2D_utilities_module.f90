!cell_mesh2d utilities module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 0.6
!Updated 12-02-2024

!Module
module cellmesh2d_utilities_mod
use cellmesh2d_data_mod
use ieee_arithmetic 
contains 


!Matrix Inverse Suborutine ===================================================
subroutine matinv(a,c,n) !Based on: Inverse matrix by Alex G. (December 2009) (a = input matrix | c = inv(a) | n = dimension of a | a is destroyed)
implicit none 

!Variables - Import
integer(in) :: n
real(dp) :: a(n,n), c(n,n)

!Variables - Local 
real(dp) :: coeff
integer(in) :: i,j,k
real(dp) :: b(n),d(n),x(n)
real(dp) :: L(n,n),U(n,n)

! step 0: initialization for matrices L and U and b
L(:,:) = 0.0d0
U(:,:) = 0.0d0
b(:) = 0.0d0

! step 1: forward elimination
do k=1,n-1
    do i=k+1,n
        coeff = a(i,k)/a(k,k)
        L(i,k) = coeff
        do j=k+1,n
            a(i,j) = a(i,j) - coeff*a(k,j)
        end do
    end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
    L(i,i) = 1.0d0
end do

! U matrix is the upper triangular part of A
do j=1,n
    do i=1,j
    U(i,j) = a(i,j)
    end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
    b(k) = 1.0d0
    d(1) = b(1)

    ! Step 3a: Solve Ld=b using the forward substitution
    do i=2,n
        d(i) = b(i)
        do j=1,i-1
            d(i) = d(i) - L(i,j)*d(j)
        end do
    end do

    ! Step 3b: Solve Ux=d using the back substitution
    x(n)=d(n)/U(n,n)
    do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
            x(i) = x(i) - U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
    end do
    ! Step 3c: fill the solutions x(n) into column k of C
    do i=1,n
        c(i,k) = x(i)
    end do
    b(k) = 0.0d0
end do
end subroutine matinv




!Subroutine to build RBF influence matrix ===========================
subroutine build_RBF_influence(R,R_sup,Npoint,point_list,vertices,cm2dopt)
implicit none 

!Variables - Import
integer(in) :: Npoint 
integer(in), dimension(:) :: point_list
real(dp) :: R_sup
real(dp), dimension(:,:) :: R,vertices 
type(cm2d_options) :: cm2dopt

!Variables - Local
integer(in) :: ii,jj 
integer(in) :: v1,v2
real(dp) :: R_relax,maxdist
real(dp) :: mindist(Npoint)

!Populate distance matrix
R(:,:) = 0.0d0 
mindist(:) = ieee_value(1.0d0,IEEE_POSITIVE_INF)
do ii=1,Npoint
    v1 = point_list(ii)
    do jj=1,Npoint
        v2 = point_list(jj)
        R(ii,jj) = norm2(vertices(v2,:) - vertices(v1,:)) 
        if ((R(ii,jj) .GT. 0.0d0) .AND. (R(ii,jj) .LT. mindist(ii))) then 
            mindist(ii) = R(ii,jj)
        end if 
    end do 
end do 
maxdist = maxval(R(1:Npoint,1:Npoint))

!Set support radius 
R_sup = cm2dopt%RBF_rsup*maxdist

!Set relaxation radius  
R_relax = cm2dopt%RBF_relaxD*maxdist

!Evaluate RBF values for each distance entry and apply relaxation where neccesary if two points are within R_relax of each other 
do ii=1,Npoint
    do jj=1,Npoint
        if (ii .NE. jj) then !Relax off diagonals if points are two close with a quadratic penalty 
            if (R(ii,jj) .LE. R_relax) then !Apply relaxation if distance is too small and between two seperate points 
                R(ii,jj) = R(ii,jj) + cm2dopt%RBF_relaxP*((R(ii,jj) - R_relax)**2)
            end if 
        end if 
        R(ii,jj) = wendlandc2(R(ii,jj),R_sup) !Evaluate RBF influence 
    end do
end do 
return 
end subroutine build_RBF_influence



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

end module cellmesh2d_utilities_mod