!cell_mesh2d io module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 5.7
!Updated 18-09-2023

!Module
module cellmesh2d_io_mod
use cellmesh2d_data_mod
contains

!Read command arguments subroutine ===========================
subroutine get_process_arguments(cm2dopt)
implicit none

!Variables - Import
type(cm2d_options) :: cm2dopt
    
!Variables - Local 
integer(in) :: nargs
integer(in32) :: arglen,argstat

!Check and process supplied command arguments
nargs = command_argument_count()
if (nargs == 0) then !Use default paths and filename
    allocate(character(len=3) :: cm2dopt%iopath)
    cm2dopt%iopath = 'io/'
    allocate(character(len=3) :: cm2dopt%optpath)
    cm2dopt%optpath = 'io/'
    allocate(character(len=23) :: cm2dopt%surfacename)
    cm2dopt%surfacename = 'cell_mesh2d_surface.dat'
elseif (nargs == 1) then !Use default paths and specified filename
    allocate(character(len=3) :: cm2dopt%iopath)
    cm2dopt%iopath = 'io/'
    allocate(character(len=3) :: cm2dopt%optpath)
    cm2dopt%optpath = 'io/'
    call get_command_argument(number=1, length=arglen)
    allocate(character(len=arglen) :: cm2dopt%surfacename)
    call get_command_argument(number=1, value=cm2dopt%surfacename, status=argstat)
else !Use specified paths and filename 
    call get_command_argument(number=1, length=arglen)
    allocate(character(len=arglen) :: cm2dopt%surfacename)
    call get_command_argument(number=1, value=cm2dopt%surfacename, status=argstat)
    call get_command_argument(number=2, length=arglen)
    allocate(character(len=arglen) :: cm2dopt%optpath)
    call get_command_argument(number=2, value=cm2dopt%optpath, status=argstat)
    call get_command_argument(number=3, length=arglen)
    allocate(character(len=arglen) :: cm2dopt%iopath)
    call get_command_argument(number=3, value=cm2dopt%iopath, status=argstat)
end if 
return 
end subroutine get_process_arguments




!Options import subroutine ===========================
subroutine cm2d_import_options(cm2dopt)
implicit none

!Variables - Import
type(cm2d_options) :: cm2dopt

!Variables - Local 
character(len=100) :: rtemp 

!Open file
open(11,file=cm2dopt%optpath//'cell_mesh2d_options.dat')

!Import options  
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%dispt
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip

read(11,*) cm2dopt%meshtype
read(11,*) !skip
read(11,*) !skip
read(11,*) rtemp
allocate(character(len=len_trim(rtemp)) :: cm2dopt%meshinout)
cm2dopt%meshinout = rtemp(1:len_trim(rtemp))
read(11,*) !skip
read(11,*) !skip
read(11,*) rtemp
allocate(character(len=len_trim(rtemp)) :: cm2dopt%surface_dir)
cm2dopt%surface_dir = rtemp(1:len_trim(rtemp))
read(11,*) !skip
read(11,*) !skip
read(11,*) rtemp
allocate(character(len=len_trim(rtemp)) :: cm2dopt%boundary_dir)
cm2dopt%boundary_dir = rtemp(1:len_trim(rtemp))
read(11,*) !skip
read(11,*) !skip
read(11,*) rtemp
allocate(character(len=len_trim(rtemp)) :: cm2dopt%meshfrmat)
cm2dopt%meshfrmat = rtemp(1:len_trim(rtemp))
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip

read(11,*) cm2dopt%Nrefine
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%NrefineB
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%Ncell_max
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%Nrefine_flood_i
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%Nrefine_flood_f
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%Nrefine_flood_B
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%far_field_bound
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%om_offset_x
read(11,*) cm2dopt%om_offset_y
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip

read(11,*) cm2dopt%set_mbounds
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%mxmin,cm2dopt%mxmax,cm2dopt%mymin,cm2dopt%mymax
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip

read(11,*) cm2dopt%EminLength
read(11,*) !skip
read(11,*) !skip 
read(11,*) cm2dopt%CminVol
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip

read(11,*) cm2dopt%NintEmax
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%elenpad
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%intcointol
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip

read(11,*) cm2dopt%surface_type
read(11,*) !skip
read(11,*) !skip 
read(11,*) cm2dopt%surfRcurvM
read(11,*) !skip
read(11,*) !skip 
read(11,*) cm2dopt%NPsinterp
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip

read(11,*) cm2dopt%Nsstype
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%nlpflood
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%nlpsmooth
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip

read(11,*) cm2dopt%ADTpadding
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%ADTmax_depth

read(11,*) !skip
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%glink_con
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%glink_nnn
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%glink_nsmooth

read(11,*) !skip
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%set_customBCs
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%remFFzones
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%remNCzones
read(11,*) !skip
read(11,*) !skip
read(11,*) cm2dopt%remISzones

!Close file 
close(11)
return 
end subroutine cm2d_import_options




!Surface data import subroutine =========================== 
subroutine import_surface_geometry(surface_mesh,cm2dopt)
implicit none 

!Variables - Import
type(surface_data) :: surface_mesh
type(cm2d_options) :: cm2dopt

!Variables - Local
integer(in) :: ii

!Open file 
open(11,file=cm2dopt%iopath//cm2dopt%surfacename)

!Import item quantities
read(11,*) surface_mesh%nvtx,surface_mesh%nfcs

!Allocate
allocate(surface_mesh%vertices(surface_mesh%nvtx,2))
allocate(surface_mesh%faces(surface_mesh%nfcs,2))
allocate(surface_mesh%vtx_active(surface_mesh%nvtx))

!Import surface verticies
do ii=1,surface_mesh%nvtx
    read(11,*) surface_mesh%vertices(ii,:) !x || y
end do
surface_mesh%vtx_active(:) = 1

!Import surface faces
do ii=1,surface_mesh%nfcs
    read(11,*) surface_mesh%faces(ii,:) !segment v1 -> v2
end do

!Close file
close(11)
return 
end subroutine import_surface_geometry




!Custom boundary condition zone import subroutine ===========================
subroutine cm2d_import_customBC_zones(cm2dopt)
implicit none

!Variables - Import
type(cm2d_options) :: cm2dopt

!Variables - Local
integer(in) :: ii

!Read file
open(11,file=cm2dopt%optpath//'cell_mesh2d_bcond_zones.dat')
    read(11,*) cm2dopt%Nzone_cBC
    allocate(cm2dopt%BC_zone_bc(cm2dopt%Nzone_cBC))
    allocate(cm2dopt%BC_zone_coords(cm2dopt%Nzone_cBC,4))
    do ii=1,cm2dopt%Nzone_cBC
        read(11,*) cm2dopt%BC_zone_bc(ii),cm2dopt%BC_zone_coords(ii,:)
    end do 
close(11)
return 
end subroutine cm2d_import_customBC_zones




!Export volume mesh subroutine ===========================
subroutine export_volume_mesh(volume_mesh,cm2dopt)
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Local
integer(in) :: ii

!Display
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '--> writing mesh to file'
end if 

!Mesh 
open(11,file=cm2dopt%iopath//'grid') !mesh file
    write(11,'(I0,A,I0,A,I0)') volume_mesh%ncell,' ',volume_mesh%nedge,' ',volume_mesh%nvtx !Properties
    do ii=1,volume_mesh%nedge
        write(11,'(I0,A,I0,A,I0,A,I0)') volume_mesh%edge(ii,1),' ',volume_mesh%edge(ii,2),' ',&
                                        volume_mesh%edge(ii,3),' ',volume_mesh%edge(ii,4) !Edges
    end do
    do ii=1,volume_mesh%nvtx
        write(11,'(I0,A,A,A,A)') ii,' ',real2F0_Xstring(volume_mesh%vertices(ii,1),16_in),&
                                    ' ',real2F0_Xstring(volume_mesh%vertices(ii,2),16_in) !Vertices
    end do
close(11)

! !Surface link data 
! open(11,file=cm2dopt%iopath//'grid_surf_link') !mesh surface to geometry surface links  
!     do ii=1,volume_mesh%nvtx_surf
!         write(11,'(I10,A,I10,A,F18.12)') volume_mesh%surf_vtx(ii),' ',volume_mesh%surf_vtx_seg(ii),' ',&
!                                          volume_mesh%surf_vtx_segfrac(ii) !vertex | link | fraction
!     end do
! close(11)

!Display
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '    {complete}'
end if
return 
end subroutine export_volume_mesh




!SU2 mesh format export subroutine =========================
subroutine export_volume_mesh_SU2(volume_mesh,cm2dopt)
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Local
integer(in) :: ii,bb,vv
integer(in) :: NBCtype,maxbc
integer(in), dimension(:), allocatable :: bcactive,nebct

!Display
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '--> writing mesh to file'
end if 

!Open file 
open(11,file=cm2dopt%iopath//'grid.su2')

!Write dimension header 
write(11,'(A)') 'NDIME= 2'

!Write vertices 
write(11,'(A,I0)') 'NPOIN= ',volume_mesh%nvtx
do ii=1,volume_mesh%nvtx
    write(11,'(A,A,A,A)') real2F0_Xstring(volume_mesh%vertices(ii,1),16_in),' ',real2F0_Xstring(volume_mesh%vertices(ii,2),16_in) !Vertices
end do

!Write elements 
write(11,'(A,I0)') 'NELEM= ',volume_mesh%ncell
do ii=1,volume_mesh%ncell
    if (volume_mesh%cells(ii)%nvtx == 3) then !triangle element
        write(11,'(I0)',advance='no') 5
    elseif (volume_mesh%cells(ii)%nvtx == 4) then !quadrilateral element
        write(11,'(I0)',advance='no') 9
    end if 
    do vv=1,volume_mesh%cells(ii)%nvtx
        write(11,'(A,I0)',advance='no') ' ',volume_mesh%cells(ii)%vertices(vv) - 1
    end do 
    write(11,'(A)') ''
end do 

!Identify all active types of boundary condition 
NBCtype = 0 
maxbc = abs(min(minval(volume_mesh%edge(:,3)),minval(volume_mesh%edge(:,4))))
allocate(bcactive(maxbc))
bcactive(:) = 0 
do ii=1,volume_mesh%nedge
    if (volume_mesh%edge(ii,3) .LT. 0) then 
        bcactive(abs(volume_mesh%edge(ii,3))) = 1
    end if 
    if (volume_mesh%edge(ii,4) .LT. 0) then 
        bcactive(abs(volume_mesh%edge(ii,4))) = 1
    end if 
end do 
NBCtype = sum(bcactive)

!Find the number of edges in each active boundary condition 
allocate(nebct(maxbc))
nebct(:) = 0 
do bb=1,maxbc
    do ii=1,volume_mesh%nedge
        if (volume_mesh%edge(ii,3) == -bb) then 
            nebct(bb) = nebct(bb) + 1
        end if 
        if (volume_mesh%edge(ii,4) == -bb) then 
            nebct(bb) = nebct(bb) + 1
        end if
    end do 
end do 

!Write boundary condition markers 
write(11,'(A,I0)') 'NMARK= ',NBCtype
do bb=1,maxbc
    if (bcactive(bb) == 1) then     

        !Write tag of this boundary condition 
        write(11,'(A,I0)') 'MARKER_TAG= ',-bb

        !Write number of edges for this boundary condition 
        write(11,'(A,I0)') 'MARKER_ELEMS= ',nebct(bb)

        !Write edges on this boundary condition 
        do ii=1,volume_mesh%nedge
            if (volume_mesh%edge(ii,3) == -bb) then 
                write(11,'(I0,A,I0,A,I0)') 3,' ',volume_mesh%edge(ii,1) - 1,' ',volume_mesh%edge(ii,2) - 1
            end if 
            if (volume_mesh%edge(ii,4) == -bb) then 
                write(11,'(I0,A,I0,A,I0)') 3,' ',volume_mesh%edge(ii,1) - 1,' ',volume_mesh%edge(ii,2) - 1
            end if
        end do 
    end if 
end do 

!Close file
close(11)

return 
end subroutine export_volume_mesh_SU2




!Write status subroutine ===========================
subroutine export_status(cm2dopt)
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt

!Write
open(11,file=cm2dopt%iopath//'cm2d_status') 
    write(11,'(i2)') cm2dopt%cm2dfailure
close(11)
return 
end subroutine export_status




!Export PPoU volume-surface to surface interpolation structure ===========================
subroutine export_vs2s_interpstruc(volume_mesh,surface_mesh,cm2dopt)
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt
type(vol_mesh_data) :: volume_mesh
type(surface_data) :: surface_mesh

!Variables - Local 
integer(in) :: vv,ii,jj

!Display
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '--> writing volume-surface to surface interpolation structure to file'
end if 

!Open 
open(11,file=cm2dopt%iopath//'vs2s_interp') 

!Write data
write(11,'(I0,A,I0,A,I0)') surface_mesh%nvtx,' ',cm2dopt%glink_nnn,' ',2*cm2dopt%glink_nsmooth+1 !number of surface vertices | number of interpolation points per surface vertex | number of smoothing points per surface vertex
write(11,'(A)') ' '

!Write each interpolation structure for each surface vertex
do vv=1,surface_mesh%nvtx
    write(11,'(I0)') vv !surface point index 
    write(11,'(I0)') volume_mesh%vsinterp(vv)%npnts_vi !number of unique vertices in the interpolation structure
    do ii=1,2*cm2dopt%glink_nsmooth+1 !write smoothing point list for this surface point 
        write(11,'(I0,A)',advance='no') volume_mesh%vsinterp(vv)%surf_smooth_pnts(ii),' '
    end do 
    write(11,'(A)') !skip to next line 
    do ii=1,2*cm2dopt%glink_nsmooth+1 !write RBF dependance for this surface point on each surface smoothing point
        write(11,'(A,A)',advance='no') real2F0_Xstring(volume_mesh%vsinterp(vv)%surf_smoothRBF(ii),12_in),' '
    end do 
    write(11,'(A)') !skip to next line 
    do ii=1,cm2dopt%glink_nnn !write volume point list for this surface point 
        write(11,'(I0,A)',advance='no') volume_mesh%vsinterp(vv)%vol_pnts(ii),' '
    end do 
    write(11,'(A)') !skip to next line 
    do ii=1,cm2dopt%glink_nnn !write RBF dependance for this surface point on each volume point 
        write(11,'(A,A)',advance='no') real2F0_Xstring(volume_mesh%vsinterp(vv)%surf2volRBF(ii),12_in),' '
    end do 
    write(11,'(A)') !skip to next line 
    do ii=1,cm2dopt%glink_nnn !write interpolation matrix 
        do jj=1,cm2dopt%glink_nnn 
            write(11,'(A,A)',advance='no') real2F0_Xstring(volume_mesh%vsinterp(vv)%Ri(ii,jj),12_in),' '
        end do 
        write(11,'(A)') !skip to next line
    end do 
    write(11,'(A)') ' '
end do 

!Close
close(11)

!Display
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '    {complete}'
end if
return 
end subroutine export_vs2s_interpstruc




!Export TECPLOT mesh file with cell data subroutine =========================  
subroutine write_cell_dataPLT(filename,volume_mesh,cell_datap)
use ieee_arithmetic
implicit none 


!Variables - Import
real(dp), dimension(:) :: cell_datap
character(*), intent(in) :: filename
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: fh,i,nperline

!Open file
fh = 11
open(fh,file=filename//'.plt',status='unknown')

!TECPLOT formatting
write(fh,'(A)',advance="no") 'VARIABLES="X" "Y" "data"'
write(fh,*) 'ZONE T="CellData"'
write(fh,'(A)',advance="no") 'VARLOCATION=([1,2]=NODAL,[3]=CELLCENTERED)'
write(fh,*) 'ZONETYPE=FEPOLYGON'
write(fh,'(A,I8)') ' Nodes=',volume_mesh%nvtx
write(fh,'(A,I8)') ' Elements=',volume_mesh%ncell
write(fh,'(A,I8)') ' Faces=',volume_mesh%nedge
write(fh,*) 'NumConnectedBoundaryFaces=0 '
write(fh,*) 'TotalNumBoundaryConnections=0 '

! These loops are because tecplot has a maximum number of characters per line
nperline = 100
write(fh,*) ( volume_mesh%vertices(i:min(i+nperline-1,volume_mesh%nvtx),1),NEW_LINE('A') , i=1,volume_mesh%nvtx,nperline )
write(fh,*) ( volume_mesh%vertices(i:min(i+nperline-1,volume_mesh%nvtx),2),NEW_LINE('A') , i=1,volume_mesh%nvtx,nperline )
write(fh,*) ( cell_datap(i:min(i+nperline-1,volume_mesh%ncell)),NEW_LINE('A') , i=1,volume_mesh%ncell,nperline )
do i=1,volume_mesh%nedge
    write(fh,*) volume_mesh%edge(i,1),volume_mesh%edge(i,2)  
end do
write(fh,*) ( max(0,volume_mesh%edge(i:min(i+nperline-1,volume_mesh%nedge),3)),NEW_LINE('A') , i=1,volume_mesh%nedge,nperline )
write(fh,*) ( max(0,volume_mesh%edge(i:min(i+nperline-1,volume_mesh%nedge),4)),NEW_LINE('A') , i=1,volume_mesh%nedge,nperline )

!Close file
close(fh)
return 
end subroutine write_cell_dataPLT




!F0.X format with leading zero function =========================
function real2F0_Xstring(val,X) result(str)

!Result 
character(len=:), allocatable :: str

!Variables - Import 
character(len=10) :: frmtI
character(len=:), allocatable :: frmt,str_I
integer(in) :: X,len_frmt,len_str
real(dp) :: val

!Set format descriptor
write(frmtI,'(I0)') X
len_frmt = len_trim(frmtI)
allocate(character(len=len_frmt) :: frmt)
frmt = frmtI(1:len_frmt)
frmt = frmt//')'
frmt = '(F0.'//frmt

!Allocate initial character
allocate(character(len=4*X) :: str_I)

!Write data to return charachter
write(str_I,frmt) val

!Allocate return character
len_str = len_trim(str_I)
allocate(character(len=len_str) :: str)
str = str_I(1:len_str)

!Assign leading zero if required
if (str(1:1) == '.') then 
    str = '0'//str
elseif (str(1:2) == '-.') then 
    str = '-0.'//str(3:len_str)
end if 
end function real2F0_Xstring


end module cellmesh2d_io_mod