!Cell Mesh 2D Surface Mesh Module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 0.2
!Updated 21-07-2023

!Module
module cellmesh2d_surface_mod
use cellmesh2d_geometry_mod
contains 


!Orient surface mesh subroutine ===========================
subroutine orient_surface(surface_mesh,cm2dopt)
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt
type(surface_data) :: surface_mesh

!Variables - Local 
integer(in) :: ff
integer(in) :: ftemp(2)
real(dp) :: AEdge,vol_total 

!Evaluate volume of objects 
vol_total = 0.0d0 
do ff=1,surface_mesh%nfcs
    AEdge = Asegment(surface_mesh%vertices(surface_mesh%faces(ff,1),:),&
                        surface_mesh%vertices(surface_mesh%faces(ff,2),:))
    vol_total = vol_total - AEdge
end do 

!Flip if negative volume 
if (vol_total .LT. 0.0d0) then 
    do ff=1,surface_mesh%nfcs
        ftemp(:) = surface_mesh%faces(ff,:)
        surface_mesh%faces(ff,1) = ftemp(2)
        surface_mesh%faces(ff,2) = ftemp(1)
    end do 
    if (cm2dopt%dispt == 1) then
        write(*,'(A)') '    {surface mesh re-oriented for positive volume}'
    end if
end if 
return 
end subroutine orient_surface




!Evaluate surface radius of curvature subroutine =========================== 
subroutine evaluate_surf_rcurv(surface_mesh,cm2dopt)
use ieee_arithmetic
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt
type(surface_data) :: surface_mesh

!Variables - Local 
integer(in) :: vv,aa,ff 
integer(in) :: cpnt,Ninterp,vb,va
integer(in) :: interp_stencil(2*cm2dopt%NPsinterp+1)

!Allocate
allocate(surface_mesh%vtx_rcurv(surface_mesh%nvtx))
surface_mesh%vtx_rcurv(:) = 0.0d0 
allocate(surface_mesh%face_rcurv(surface_mesh%nfcs))
surface_mesh%face_rcurv(:) = 0.0d0 

!Set centre point 
cpnt = cm2dopt%NPsinterp+1

!Set Ninterp
Ninterp = 2*cm2dopt%NPsinterp+1

!Evaluate for each vertex
do vv=1,surface_mesh%nvtx

    !Build interpolation stencil  
    interp_stencil(:) = 0 
    vb = vv 
    do aa=1,cm2dopt%NPsinterp !previous vertices
        va = get_previous_vertex(surface_mesh,vb)
        interp_stencil(cpnt-aa) = va
        vb = va 
    end do 
    interp_stencil(cpnt) = vv 
    vb = vv 
    do aa=1,cm2dopt%NPsinterp !next vertices
        va = get_next_vertex(surface_mesh,vb)
        interp_stencil(cpnt+aa) = va
        vb = va 
    end do 

    !Evaluate local radius of curvature at this vertex
    surface_mesh%vtx_rcurv(vv) = surface_rcurv(Ninterp,interp_stencil,surface_mesh%vertices)
    surface_mesh%vtx_rcurv(vv) = surface_mesh%vtx_rcurv(vv)/cm2dopt%surfRcurvM
end do 

!Project to faces
do ff=1,surface_mesh%nfcs
    surface_mesh%face_rcurv(ff) = 0.5d0*(surface_mesh%vtx_rcurv(surface_mesh%faces(ff,1)) + &
    surface_mesh%vtx_rcurv(surface_mesh%faces(ff,2)))
end do 

!Debug 
! open(11,file='io/rcurv.dat')
! do vv=1,surface_mesh%nvtx
!     write(11,*) surface_mesh%vtx_rcurv(vv)
! end do 
! close(11)
return 
end subroutine evaluate_surf_rcurv




!V2F construction subroutine ===========================
subroutine construct_surface_v2f(surface_mesh,cm2dopt)
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt
type(surface_data) :: surface_mesh

!Variables - Local 
integer(in) :: ii 
integer(in) :: vtx_nseg(surface_mesh%nvtx)

!Construct 
vtx_nseg(:) = 0 
allocate(surface_mesh%v2f(surface_mesh%nvtx,2))
surface_mesh%v2f(:,:) = 0 
do ii=1,surface_mesh%nfcs
    vtx_nseg(surface_mesh%faces(ii,1)) = vtx_nseg(surface_mesh%faces(ii,1)) + 1
    vtx_nseg(surface_mesh%faces(ii,2)) = vtx_nseg(surface_mesh%faces(ii,2)) + 1
    if ((vtx_nseg(surface_mesh%faces(ii,1)) .GT. 2) .OR. (vtx_nseg(surface_mesh%faces(ii,2)) .GT. 2)) then 
        cm2dopt%cm2dfailure = 1
        print *, '** non-manifold input surface geometry detected at segment : ', ii
        return 
    end if
    surface_mesh%v2f(surface_mesh%faces(ii,1),2) = ii
    surface_mesh%v2f(surface_mesh%faces(ii,2),1) = ii
end do 
return 
end subroutine construct_surface_v2f




!Get next vertex function ===========================
function get_next_vertex(surface_mesh,vb) result(vn)
implicit none 

!Variables - Import
integer(in) :: vb,vn
type(surface_data) :: surface_mesh

!Variables - Local 
integer(in) :: fn 

!Get next vertex
fn = surface_mesh%v2f(vb,2)
if (surface_mesh%faces(fn,1) == vb) then 
    vn = surface_mesh%faces(fn,2)
else
    vn = surface_mesh%faces(fn,1)
end if
return 
end function get_next_vertex




!Get previous vertex function ===========================
function get_previous_vertex(surface_mesh,vb) result(vp)
implicit none 

!Variables - Import
integer(in) :: vb,vp
type(surface_data) :: surface_mesh

!Variables - Local 
integer(in) :: fp 

!Get previous vertex
fp = surface_mesh%v2f(vb,1)
if (surface_mesh%faces(fp,1) == vb) then 
    vp = surface_mesh%faces(fp,2)
else
    vp = surface_mesh%faces(fp,1)
end if
return 
end function get_previous_vertex


end module cellmesh2d_surface_mod