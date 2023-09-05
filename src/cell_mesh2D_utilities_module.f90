!cell_mesh2d utilities module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 0.51
!Updated 19-07-2023

!Module
module cellmesh2d_utilities_mod
use cellmesh2d_data_mod
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


end module cellmesh2d_utilities_mod