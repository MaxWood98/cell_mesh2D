!cell_mesh2d minD mesh module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 0.1
!Updated 16-02-2023

!Module 
module cellmesh2d_mind_mod
use cellmesh2d_distF_mod
contains 


!MinD mesh construction subroutine ================================================
subroutine build_minD_Omesh(volume_mesh,surface_mesh,cm2dopt)
implicit none 

!Variables - Import 
type(cm2d_options) :: cm2dopt
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ee,vv
integer(in) :: BCactive(volume_mesh%nedge)
real(dp) :: cell_d(volume_mesh%ncell),vtx_d(volume_mesh%nvtx),grad2_d(volume_mesh%ncell)
real(dp) :: grad_d(volume_mesh%ncell,2)
type(meshdata) :: meshdat

!Preprocess mesh 
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '--> preprocessing mesh for distance field construction'
end if
call preprocess_mesh(meshdat,volume_mesh)

!Construct distance field 
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '--> constructing distance field'
end if
BCactive(:) = 0
do ee=1,volume_mesh%nedge !set wall boundary conditions to active 
    if ((volume_mesh%edge(ee,3) == -1) .OR. (volume_mesh%edge(ee,4) == -1)) then 
        BCactive(ee) = 1
    elseif ((volume_mesh%edge(ee,3) == -2) .OR. (volume_mesh%edge(ee,4) == -2)) then 
        ! BCactive(ee) = 1 
    end if 
end do 
call build_distance_field(cell_d,grad_d,grad2_d,volume_mesh,meshdat,BCactive,cm2dopt)
vtx_d = project_cell2vtx(cell_d,meshdat,volume_mesh,1_in) 

!Build mesh 





open(11,file='io/vtxd')
do vv=1,volume_mesh%nvtx
    write(11,*) vtx_d(vv)
end do 
close(11)


! open(11,file='io/celld')
! do vv=1,volume_mesh%ncell
!     write(11,*) cell_d(vv)
! end do 
! close(11)

call saveVolPLT('io/cell_dist.plt',volume_mesh,cell_d,grad_d,&
(grad_d(:,1) + grad_d(:,2))*(grad_d(:,1) - grad_d(:,2)))

call saveVolPLT('io/cell_dist_g2.plt',volume_mesh,cell_d,grad_d,abs(grad2_d))



! call saveVolPLT('io/cell_dist.plt',volume_mesh,cell_d,grad_d,&
! abs((grad_d(:,1) - grad_d(:,2))*(grad_d(:,2) - grad_d(:,1))))


! call saveVolPLT('io/cell_dist.plt',volume_mesh,cell_d,grad_d,&
!  abs((grad_d(:,1) - grad_d(:,2)))*abs((grad_d(:,2) - grad_d(:,1)))*(grad_d(:,2) + grad_d(:,1))   )

! call saveVolPLT('io/cell_dist.plt',volume_mesh,cell_d,grad_d,&
! (grad_d(:,2) + grad_d(:,1))*(grad_d(:,1) - grad_d(:,2))  )


! call saveVolPLT('io/cell_dist.plt',volume_mesh,cell_d,grad_d,abs(grad2_d) )

return 
end subroutine build_minD_Omesh






!Export TECPLOT solution file subroutine ========================= (from edge2d by Tom Rendall)
subroutine saveVolPLT(filename,volume_mesh,cell_dist,grad_dist,propT)
implicit none 

!Variables - Import
character(*), intent(in) :: filename
real(dp), dimension(:) :: cell_dist,propT
real(dp), dimension(:,:) :: grad_dist
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: fh,i,nperline
real(dp) :: vtx_x(volume_mesh%nvtx),vtx_y(volume_mesh%nvtx)

!Open file
fh = 11
open(fh,file=filename,status='unknown')

!TECPLOT formatting
write(fh,'(A)',advance="no") 'VARIABLES="X" "Y"'
write(fh,*) '"distance" "gradx" "grady" "property"'
write(fh,*) 'ZONE T="DistanceField"'
write(fh,'(A)',advance="no") 'VARLOCATION=([1,2]=NODAL'
write(fh,*) ',[3,4,5,6]=CELLCENTERED)'
write(fh,*) 'ZONETYPE=FEPOLYGON'
write(fh,'(A,I8)') ' Nodes=',volume_mesh%nvtx
write(fh,'(A,I8)') ' Elements=',volume_mesh%ncell
write(fh,'(A,I8)') ' Faces=',volume_mesh%nedge
write(fh,*) 'NumConnectedBoundaryFaces=0 '
write(fh,*) 'TotalNumBoundaryConnections=0 '

! These loops are because tecplot has a maximum number of characters per line
nperline = 100
vtx_x(:) = volume_mesh%vertices(:,1)
vtx_y(:) = volume_mesh%vertices(:,2)
write(fh,*) ( vtx_x(i:min(i+nperline-1,volume_mesh%nvtx)),NEW_LINE('A') , i=1,volume_mesh%nvtx,nperline )
write(fh,*) ( vtx_y(i:min(i+nperline-1,volume_mesh%nvtx)),NEW_LINE('A') , i=1,volume_mesh%nvtx,nperline )
write(fh,*) ( cell_dist(i:min(i+nperline-1,volume_mesh%ncell)),NEW_LINE('A') , i=1,volume_mesh%ncell,nperline )
write(fh,*) ( grad_dist(i:min(i+nperline-1,volume_mesh%ncell),1),NEW_LINE('A') , i=1,volume_mesh%ncell,nperline )
write(fh,*) ( grad_dist(i:min(i+nperline-1,volume_mesh%ncell),2),NEW_LINE('A') , i=1,volume_mesh%ncell,nperline )
write(fh,*) ( propT(i:min(i+nperline-1,volume_mesh%ncell)),NEW_LINE('A') , i=1,volume_mesh%ncell,nperline )
do i=1,volume_mesh%nedge
    write(fh,*) volume_mesh%edge(i,1),volume_mesh%edge(i,2)  
end do
write(fh,*) ( max(0,volume_mesh%edge(i:min(i+nperline-1,volume_mesh%nedge),3)),NEW_LINE('A') , i=1,volume_mesh%nedge,nperline )
write(fh,*) ( max(0,volume_mesh%edge(i:min(i+nperline-1,volume_mesh%nedge),4)),NEW_LINE('A') , i=1,volume_mesh%nedge,nperline )


!Close file
close(fh)
return 
end subroutine saveVolPLT


end module cellmesh2d_mind_mod