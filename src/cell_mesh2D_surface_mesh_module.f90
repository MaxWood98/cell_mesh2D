!Cell Mesh 2D Surface Mesh Module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 0.5
!Updated 07-11-2023

!Module
module cellmesh2d_surface_mod
use cellmesh2d_geometry_mod
contains 


!Preprocess surface mesh subroutine ===========================
subroutine preprocess_surface_mesh(surface_mesh,cm2dopt)
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt
type(surface_data) :: surface_mesh

!Variables - Local 
integer(in) :: ii
integer(in) :: nrem

!Clean surface mesh to remove zero length segments 
do ii=1,10
    call clean_surface(surface_mesh,cm2dopt,nrem)
    if (nrem == 0) then 
        exit 
    end if 
end do 

!Orient surface mesh for positive object volume (this ensures normal vector convention is correct)
call orient_surface(surface_mesh,cm2dopt)

!Build v2f for surface geometry and check for non-manifold surface geometry 
call construct_surface_v2f(surface_mesh,cm2dopt)
if (cm2dopt%cm2dfailure == 1) then
    return 
end if 

!Evaluate surface mesh segment curvature 
call evaluate_surf_rcurv(surface_mesh,cm2dopt)
return 
end subroutine preprocess_surface_mesh




!Clean surface with retained vertex indexing subroutine ===========================
subroutine clean_surface(surface_mesh,cm2dopt,nrem)
implicit none 

!Variables - Import
integer(in) :: nrem
type(cm2d_options) :: cm2dopt
type(surface_data) :: surface_mesh

!Variables - Local 
integer(in) :: ii,v1,v2,NfaceN,NfaceRem
integer(in) :: faceindexN(surface_mesh%nfcs)
integer(in) :: vtx_idxN(surface_mesh%nvtx)
integer(in) :: faces_temp(surface_mesh%nfcs,2)
real(dp) :: lenseg(surface_mesh%nfcs)

!Initialise 
nrem = 0 

!Evaluate segment norms 
lenseg(:) = 0.0d0 
do ii=1,surface_mesh%nfcs
    v1 = surface_mesh%faces(ii,1)
    v2 = surface_mesh%faces(ii,2)
    lenseg(ii) = norm2(surface_mesh%vertices(v1,:) - surface_mesh%vertices(v2,:))
end do 

!Index faces to retain 
NfaceN = 0 
NfaceRem = 0 
faceindexN(:) = 0 
do ii=1,surface_mesh%nfcs
    if (lenseg(ii) .GT. 0.0d0) then 
        NfaceN = NfaceN + 1
        faceindexN(ii) = NfaceN
    else 
        NfaceRem = NfaceRem + 1
    endif 
end do 

!Return if no faces need to be removed 
if (NfaceRem == 0) then 
    return 
end if 

!Merge vertex indecies in faces where they have been merged 
vtx_idxN(:) = 0 
do ii=1,surface_mesh%nfcs
    v1 = surface_mesh%faces(ii,1)
    v2 = surface_mesh%faces(ii,2)
    if (faceindexN(ii) == 0) then !removed face -> merge vertices
        if ((vtx_idxN(v1) == 0) .AND. (vtx_idxN(v2) == 0)) then 
            vtx_idxN(v1) = v1
            vtx_idxN(v2) = vtx_idxN(v1)
        elseif ((vtx_idxN(v1) .NE. 0) .AND. (vtx_idxN(v2) == 0)) then 
            vtx_idxN(v2) = vtx_idxN(v1)
        elseif ((vtx_idxN(v1) == 0) .AND. (vtx_idxN(v2) .NE. 0)) then 
            vtx_idxN(v1) = vtx_idxN(v2)
        else
            ! print *, 'ecase'
        end if 
    end if
end do 
do ii=1,surface_mesh%nvtx
    if (vtx_idxN(ii) == 0) then 
        vtx_idxN(ii) = ii 
    end if 
end do 

!Build new faces
faces_temp(:,:) = surface_mesh%faces(:,:)
deallocate(surface_mesh%faces)
allocate(surface_mesh%faces(NfaceN,2))
do ii=1,surface_mesh%nfcs
    if (faceindexN(ii) .GT. 0) then 
        surface_mesh%faces(faceindexN(ii),:) = faces_temp(ii,:)
    end if 
end do 
surface_mesh%nfcs = NfaceN

!Remap face vertex indecies 
do ii=1,surface_mesh%nfcs
    surface_mesh%faces(ii,1) = vtx_idxN(surface_mesh%faces(ii,1))
    surface_mesh%faces(ii,2) = vtx_idxN(surface_mesh%faces(ii,2))
    ! print *, surface_mesh%faces(ii,:)
end do 

!Set active vertices 
surface_mesh%vtx_active(:) = 0 
do ii=1,surface_mesh%nfcs
    surface_mesh%vtx_active(surface_mesh%faces(ii,:)) = 1
end do 

!Store number of removed faces
nrem = NfaceRem

!Display
if (cm2dopt%dispt == 1) then
    write(*,'(A,I0,A)') '    {removed ',NfaceRem,' zero length faces}'
end if
return 
end subroutine clean_surface




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
    if (surface_mesh%vtx_active(vv) == 1) then 

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
    end if 
end do 

!Project to faces
do ff=1,surface_mesh%nfcs
    surface_mesh%face_rcurv(ff) = 0.5d0*(surface_mesh%vtx_rcurv(surface_mesh%faces(ff,1)) + &
    surface_mesh%vtx_rcurv(surface_mesh%faces(ff,2)))
    ! print *, surface_mesh%face_rcurv(ff)
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
! do ii=1,surface_mesh%nvtx
!     print *, surface_mesh%v2f(ii,:)
! end do 
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