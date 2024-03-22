!cell_mesh2d io module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 7.2
!Updated 22-03-2024

!Module
module cellmesh2d_io_mod
use io_utilities
contains

!Read command arguments subroutine ===========================
subroutine get_process_arguments(cm2dopt)
implicit none

!Variables - Import
type(cm2d_options) :: cm2dopt

!Variables - Local 
integer(in32) :: ii,jj
integer(in32) :: nargs,pathpos
character(len=:), allocatable :: argcurr,pathtemp


!Check and process supplied command arguments 
nargs = command_argument_count()
if (nargs == 0) then 
    write(*,'(A)') '** at least one argument [mode] must be supplied,& 
    & optionaly followed by the paths of input and output items'
    print *, '** paths: '
    print *, '-o [options file with path]'
    print *, '-s [surface file with path]'
    stop
else

    !Get operation mode 
    cm2dopt%mode = get_command_argument_n_str(1)

    !Error if incorrect mode specified 
    if (cm2dopt%mode == 'check') then 
        !do nothing
    elseif (cm2dopt%mode == 'mesh') then 
        !do nothing
    elseif (cm2dopt%mode == 'project') then 
        !do nothing
    else
        write(*,'(A,A)') '** unknown mode option requested : ',cm2dopt%mode
        stop
    end if 

     !Scan for additional arguments 
    do ii=2,nargs-1

        !Current argument 
        argcurr = get_command_argument_n_str(ii)

        !If option tag
        if (argcurr == '-o') then !next is options filename with path 

            !Read options filepath and extract optpath
            cm2dopt%optpath = get_command_argument_n_str(ii+1)
            pathpos = 0 
            do jj=len_trim(cm2dopt%optpath),1,-1
                if ((cm2dopt%optpath(jj:jj) == '/') .OR. (cm2dopt%optpath(jj:jj) == '\')) then 
                    pathpos = jj 
                    exit 
                end if 
            end do 
            if (pathpos .NE. 0) then !path
                pathtemp = cm2dopt%optpath
                cm2dopt%options_filename = pathtemp(pathpos+1:len_trim(pathtemp))
                cm2dopt%optpath = pathtemp(1:pathpos)
            else !no path 
                cm2dopt%options_filename = cm2dopt%optpath
                cm2dopt%optpath = ''
            end if 
            ! print *, pathpos
            ! print *, cm2dopt%optpath
            ! print *, cm2dopt%options_filename
        elseif (argcurr == '-s') then !next is surface filename with path

            !Read surface filepath and extract iopath 
            cm2dopt%iopath = get_command_argument_n_str(ii+1)
            pathpos = 0 
            do jj=len_trim(cm2dopt%iopath),1,-1
                if ((cm2dopt%iopath(jj:jj) == '/') .OR. (cm2dopt%iopath(jj:jj) == '\')) then 
                    pathpos = jj 
                    exit 
                end if 
            end do 
            if (pathpos .NE. 0) then !path
                pathtemp = cm2dopt%iopath
                cm2dopt%surface_filename = pathtemp(pathpos+1:len_trim(pathtemp))
                cm2dopt%iopath = pathtemp(1:pathpos)
            else !no path 
                cm2dopt%surface_filename = cm2dopt%iopath
                cm2dopt%iopath = ''
            end if 
            ! print *, pathpos
            ! print *, cm2dopt%iopath
            ! print *, cm2dopt%surface_filename
        end if 
    end do 
end if 
return 
end subroutine get_process_arguments




!Set default paths subroutine ===========================
subroutine set_default_paths(cm2dopt)
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt

!Set paths
cm2dopt%iopath = ''
cm2dopt%optpath = ''
cm2dopt%surface_filename = 'cell_mesh2d_surface'
cm2dopt%options_filename = 'cell_mesh2d_options'
cm2dopt%bcondzone_filename = 'cell_mesh2d_bcond_zones'
return 
end subroutine set_default_paths




!Set default options subroutine ===========================
subroutine set_default_options(cm2dopt)
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt

!Set general options 
cm2dopt%dispt = .true.
cm2dopt%meshtype = 'cutcell'
cm2dopt%meshinout = 'out'
cm2dopt%surface_dir = 'in'
cm2dopt%boundary_dir = 'in'

!Set quadtree refinement options 
cm2dopt%Nrefine = 11
cm2dopt%NrefineB = 0
cm2dopt%Ncell_max = 200000
cm2dopt%Nrefine_flood_i = 10
cm2dopt%Nrefine_flood_f = 10
cm2dopt%Nrefine_flood_B = 2
cm2dopt%far_field_bound = 15.0d0 
cm2dopt%om_offset_x = 0.0d0 
cm2dopt%om_offset_y = 0.0d0 

!Set domain bound options
cm2dopt%set_mbounds = 'no'
cm2dopt%mxmin = -10.0d0
cm2dopt%mxmax = 10.0d0
cm2dopt%mymin = -10.0d0 
cm2dopt%mymax = 10.0d0 

!Set mesh cleaning options
cm2dopt%EminLength = 1e-8
cm2dopt%sEminLength = 0.005d0 
cm2dopt%CminVol = 0.1d0 
cm2dopt%srfinclean = 'no'

!Set geometry intersection options
cm2dopt%NintEmax = 50
cm2dopt%elenpad = 0.0d0 
cm2dopt%intcointol = 1e-8

!Set mesh surface options 
cm2dopt%surface_type = 'simplified'
cm2dopt%surfRcurvM = 1.0d0 
cm2dopt%vtx_sharp_dpval = 0.5d0 
cm2dopt%NPsinterp = 20 

!Set inflation layer options
cm2dopt%build_inflayer = 'no'
cm2dopt%inflayer_height = 0.075d0 
cm2dopt%inflayer_nlayer = 0
cm2dopt%inflayer_nintstep_max = 10000
cm2dopt%inflayer_h0 = 0.01d0 
cm2dopt%inflayer_nbclinesearch = 100
cm2dopt%inflayer_dvw = 0.1d0 
cm2dopt%inflayer_cvxdp = 0.05d0
cm2dopt%inflayer_ew = 0.1d0 
cm2dopt%inflayer_ebcbase = 0.1d0 
cm2dopt%inflayer_enormw = 0.5d0 
cm2dopt%inflayer_enflood = 50
cm2dopt%inflayer_ensubiter = 10
cm2dopt%inflayer_cvxep = 1.0d0 
cm2dopt%inflayer_lreb = 1.0d0 
cm2dopt%inflayer_stepnacheck = 2

!Set mesh smoothing options 
cm2dopt%Nsstype = 'none'
cm2dopt%nlpflood = 5
cm2dopt%nlpsmooth = 10

!Set ADtree options
cm2dopt%ADTpadding = 0.0d0
cm2dopt%ADTmax_depth = 10 
cm2dopt%ADTminNodedivsize = 10

!Set gradient projection options 
cm2dopt%glink_con_exp = 'no' 
cm2dopt%glink_type = 'rbf'
cm2dopt%glink_nnn = 10
cm2dopt%glink_nsmooth = 0 

!Set radial basis function options 
cm2dopt%RBF_rsup = 50.0d0 
cm2dopt%RBF_relaxD = 0.05d0 
cm2dopt%RBF_relaxP = 0.5d0 

!Set boundary condition options 
cm2dopt%set_customBCs = 'no' 
cm2dopt%remFFzones = 'no'
cm2dopt%remNCzones = 'no' 
cm2dopt%remISzones = 'no'

!Set distance field options 
cm2dopt%dfield_niter = 10000
cm2dopt%dfield_cfl = 0.5d0 
cm2dopt%dfield_kd = 0.05d0 
cm2dopt%dfield_cres = -8.0d0 
return 
end subroutine set_default_options




!Options import subroutine ===========================
subroutine cm2d_import_options(cm2dopt)
implicit none

!Variables - Import
type(cm2d_options) :: cm2dopt

!Check if file exists
if (.NOT.file_exists(cm2dopt%optpath//cm2dopt%options_filename)) then 
    write(*,'(A)') '    ** cannot locate options file: '//cm2dopt%optpath//cm2dopt%options_filename
    stop
end if 

!Open file
open(11,file=cm2dopt%optpath//cm2dopt%options_filename)

!Set general options 
call set_log_opt(cm2dopt%dispt,11,'condisp')
call set_str_opt(cm2dopt%meshtype,11,'meshtype')
call set_str_opt(cm2dopt%meshinout,11,'meshinout')
call set_str_opt(cm2dopt%surface_dir,11,'surfnormdir')
call set_str_opt(cm2dopt%boundary_dir,11,'bndrynormdir')

!Set quadtree refinement options 
call set_int_opt(cm2dopt%Nrefine,11,'nqtrefine')
call set_int_opt(cm2dopt%NrefineB,11,'nboostqtrefine')
call set_int_opt(cm2dopt%Ncell_max,11,'ncellmax')
call set_int_opt(cm2dopt%Nrefine_flood_i,11,'nadjfloodi')
call set_int_opt(cm2dopt%Nrefine_flood_f,11,'nadjfloodf')
call set_int_opt(cm2dopt%Nrefine_flood_B,11,'nadjfloodb')
call set_real_opt(cm2dopt%far_field_bound,11,'farfielddist')
call set_real_opt(cm2dopt%om_offset_x,11,'offsett_x')
call set_real_opt(cm2dopt%om_offset_y,11,'offsett_y')

!Set domain bound options
call set_str_opt(cm2dopt%set_mbounds,11,'forcebounds')
call set_real_opt(cm2dopt%mxmin,11,'bound_xmin')
call set_real_opt(cm2dopt%mxmax,11,'bound_xmax')
call set_real_opt(cm2dopt%mymin,11,'bound_ymin')
call set_real_opt(cm2dopt%mymax,11,'bound_ymax')

!Set mesh cleaning options
call set_real_opt(cm2dopt%EminLength,11,'eminlength')
call set_real_opt(cm2dopt%sEminLength,11,'srfeminlength')
call set_real_opt(cm2dopt%CminVol,11,'cminvol')
call set_str_opt(cm2dopt%srfinclean,11,'srfinclean')

!Set geometry intersection options
call set_int_opt(cm2dopt%NintEmax,11,'enintmax')
call set_real_opt(cm2dopt%elenpad,11,'eintpad')
call set_real_opt(cm2dopt%intcointol,11,'intcointol')

!Set mesh surface options 
call set_str_opt(cm2dopt%surface_type,11,'surftype')
call set_real_opt(cm2dopt%surfRcurvM,11,'scurvmult')
call set_int_opt(cm2dopt%NPsinterp,11,'scurvnpnt')

!Inflation layer options 
call set_str_opt(cm2dopt%build_inflayer,11,'build_inflayer')
call set_real_opt(cm2dopt%inflayer_height,11,'inflayer_height')
call set_int_opt(cm2dopt%inflayer_nlayer,11,'inflayer_nlayer')
call set_real_opt(cm2dopt%inflayer_h0,11,'inflayer_h0')
call set_int_opt(cm2dopt%inflayer_nintstep_max,11,'inflayer_nintstep_max')
call set_real_opt(cm2dopt%inflayer_dvw,11,'inflayer_wd')
call set_real_opt(cm2dopt%inflayer_cvxdp,11,'inflayer_cvxdp')
call set_real_opt(cm2dopt%inflayer_ew,11,'inflayer_we')
call set_real_opt(cm2dopt%inflayer_ebcbase,11,'inflayer_ebcbase')
call set_real_opt(cm2dopt%inflayer_enormw,11,'inflayer_enormw')
call set_int_opt(cm2dopt%inflayer_enflood,11,'inflayer_enflood')
call set_int_opt(cm2dopt%inflayer_ensubiter,11,'inflayer_ensubiter')
call set_real_opt(cm2dopt%inflayer_cvxep,11,'inflayer_cvxep')
call set_real_opt(cm2dopt%inflayer_lreb,11,'inflayer_lreb')
call set_int_opt(cm2dopt%inflayer_nbclinesearch,11,'inflayer_nbclinesearch')
call set_int_opt(cm2dopt%inflayer_stepnacheck,11,'inflayer_stepnacheck')

!Set mesh smoothing options 
call set_str_opt(cm2dopt%Nsstype,11,'smoothingtype')
call set_int_opt(cm2dopt%nlpflood,11,'smoothingnflood')
call set_int_opt(cm2dopt%nlpsmooth,11,'smoothingniter')

!Set ADtree options
call set_real_opt(cm2dopt%ADTpadding,11,'adtpadding')
call set_int_opt(cm2dopt%ADTmax_depth,11,'adtndimcycle')
call set_int_opt(cm2dopt%ADTminNodedivsize,11,'adtmindivnsize')

!Set gradient projection options 
call set_str_opt(cm2dopt%glink_con_exp,11,'glinkconstructexp')
call set_str_opt(cm2dopt%glink_type,11,'glinktype')
call set_int_opt(cm2dopt%glink_nnn,11,'glinkrbfnpnt')
call set_int_opt(cm2dopt%glink_nsmooth,11,'glinknpntsmooth')

!Set radial basis function options 
call set_real_opt(cm2dopt%RBF_rsup,11,'rbfrsup')
call set_real_opt(cm2dopt%RBF_relaxD,11,'rbfrrelax')
call set_real_opt(cm2dopt%RBF_relaxP,11,'rbfrelaxp')

!Set boundary condition options 
call set_str_opt(cm2dopt%set_customBCs,11,'setcustombcs')
call set_str_opt(cm2dopt%remFFzones,11,'remffczones')
call set_str_opt(cm2dopt%remNCzones,11,'remncczones')
call set_str_opt(cm2dopt%remISzones,11,'remwallonlyzones')

!Set distance field options 


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

!Check if file exists
if (.NOT.file_exists(cm2dopt%iopath//cm2dopt%surface_filename)) then 
    write(*,'(A)') '    ** cannot locate surface file: '//cm2dopt%iopath//cm2dopt%surface_filename
    stop
end if 

!Open file 
open(11,file=cm2dopt%iopath//cm2dopt%surface_filename)

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

!Check if file exists
if (.NOT.file_exists(cm2dopt%optpath//cm2dopt%bcondzone_filename)) then 
    write(*,'(A)') '    ** cannot locate boundary contition file: '//cm2dopt%optpath//cm2dopt%bcondzone_filename
    stop
end if 

!Read file
open(11,file=cm2dopt%optpath//cm2dopt%bcondzone_filename)
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
if (cm2dopt%dispt) then
    write(*,'(A)') '--> writing mesh to file'
end if 

!Mesh 
open(11,file=cm2dopt%iopath//'grid') !mesh file
    write(11,'(I0,A,I0,A,I0)') volume_mesh%ncell,' ',volume_mesh%nedge,' ',volume_mesh%nvtx !Properties
    do ii=1,volume_mesh%nedge !Edges
        write(11,'(I0,A,I0,A,I0,A,I0)') volume_mesh%edge(ii,1),' ',volume_mesh%edge(ii,2),' ',&
                                        volume_mesh%edge(ii,3),' ',volume_mesh%edge(ii,4) !Edges
    end do
    do ii=1,volume_mesh%nvtx !Vertices
        write(11,'(I0,A,A,A,A)') ii,' ',real2F0_Xstring(volume_mesh%vertices(ii,1),16_in),&
                                    ' ',real2F0_Xstring(volume_mesh%vertices(ii,2),16_in) !Vertices
    end do
    write(11,'(I0)') volume_mesh%nvtx_surf
    do ii=1,volume_mesh%nvtx_surf !Surface link data 
        write(11,'(I0,A,I0,A,A)') volume_mesh%surf_vtx(ii),' ',volume_mesh%surf_vtx_seg(ii),' ',&
        real2F0_Xstring(volume_mesh%surf_vtx_segfrac(ii),16_in) !vertex | link | fraction
    end do
close(11)

! !Surface link data 
! open(11,file=cm2dopt%iopath//'grid_surf_link') !mesh surface to geometry surface links  
!     do ii=1,volume_mesh%nvtx_surf
!         write(11,'(I0,A,I0,A,A)') volume_mesh%surf_vtx(ii),' ',volume_mesh%surf_vtx_seg(ii),' ',&
!         real2F0_Xstring(volume_mesh%surf_vtx_segfrac(ii),16_in) !vertex | link | fraction
!     end do
! close(11)

!Display
if (cm2dopt%dispt) then
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
integer(in) :: ii,bb,vv,ee
integer(in) :: NBCtype,maxbc,v1,v2
integer(in), dimension(:), allocatable :: bcactive,nebct

!Display
if (cm2dopt%dispt) then
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
    else
        write(*,'(A,I0,A,I0)') '    ** warning: non tri/quad cell detected with ',volume_mesh%cells(ii)%nvtx,' vertices at cell ',ii
    end if 
    do vv=1,volume_mesh%cells(ii)%nvtx
        write(11,'(A,I0)',advance='no') ' ',volume_mesh%cells(ii)%vertices(vv) - 1
    end do 
    write(11,'(A)') ''
end do 

!Identify all active types of boundary condition 
NBCtype = 0 
maxbc = 0
do ii=1,volume_mesh%ncell
    do ee=1,volume_mesh%cells(ii)%nedge
        if (volume_mesh%cells(ii)%boundary_condition(ee) .LT. maxbc) then 
            maxbc = volume_mesh%cells(ii)%boundary_condition(ee)
        end if 
    end do 
end do 
maxbc = abs(maxbc)
allocate(bcactive(maxbc))
bcactive(:) = 0 
do ii=1,volume_mesh%ncell
    do ee=1,volume_mesh%cells(ii)%nedge
        if (volume_mesh%cells(ii)%boundary_condition(ee) .LT. 0) then 
            bcactive(abs(volume_mesh%cells(ii)%boundary_condition(ee))) = 1
        end if 
    end do 
end do 
NBCtype = sum(bcactive)

!Find the number of edges in each active boundary condition 
allocate(nebct(maxbc))
nebct(:) = 0 
do bb=1,maxbc
    do ii=1,volume_mesh%ncell
        do ee=1,volume_mesh%cells(ii)%nedge
            if (volume_mesh%cells(ii)%boundary_condition(ee) == -bb) then 
                nebct(bb) = nebct(bb) + 1
            end if 
        end do 
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
        do ii=1,volume_mesh%ncell
            do ee=1,volume_mesh%cells(ii)%nedge
                if (volume_mesh%cells(ii)%boundary_condition(ee) == -bb) then 

                    !Vertices on this edge
                    v1 = ee 
                    v2 = modulo(ee,volume_mesh%cells(ii)%nedge) + 1
                    v1 = volume_mesh%cells(ii)%vertices(v1)
                    v2 = volume_mesh%cells(ii)%vertices(v2)

                    !Write
                    write(11,'(I0,A,I0,A,I0)') 3,' ',v1 - 1,' ',v2 - 1
                end if 
            end do 
        end do 
    end if 
end do 















! !Identify all active types of boundary condition 
! NBCtype = 0 
! maxbc = abs(min(minval(volume_mesh%edge(:,3)),minval(volume_mesh%edge(:,4))))
! allocate(bcactive(maxbc))
! bcactive(:) = 0 
! do ii=1,volume_mesh%nedge
!     if (volume_mesh%edge(ii,3) .LT. 0) then 
!         bcactive(abs(volume_mesh%edge(ii,3))) = 1
!     end if 
!     if (volume_mesh%edge(ii,4) .LT. 0) then 
!         bcactive(abs(volume_mesh%edge(ii,4))) = 1
!     end if 
! end do 
! NBCtype = sum(bcactive)

! !Find the number of edges in each active boundary condition 
! allocate(nebct(maxbc))
! nebct(:) = 0 
! do bb=1,maxbc
!     do ii=1,volume_mesh%nedge
!         if (volume_mesh%edge(ii,3) == -bb) then 
!             nebct(bb) = nebct(bb) + 1
!         end if 
!         if (volume_mesh%edge(ii,4) == -bb) then 
!             nebct(bb) = nebct(bb) + 1
!         end if
!     end do 
! end do 

! !Write boundary condition markers 
! write(11,'(A,I0)') 'NMARK= ',NBCtype
! do bb=1,maxbc
!     if (bcactive(bb) == 1) then     

!         !Write tag of this boundary condition 
!         write(11,'(A,I0)') 'MARKER_TAG= ',-bb

!         !Write number of edges for this boundary condition 
!         write(11,'(A,I0)') 'MARKER_ELEMS= ',nebct(bb)

!         !Write edges on this boundary condition 
!         do ii=1,volume_mesh%nedge
!             if (volume_mesh%edge(ii,3) == -bb) then 
!                 write(11,'(I0,A,I0,A,I0)') 3,' ',volume_mesh%edge(ii,1) - 1,' ',volume_mesh%edge(ii,2) - 1
!             end if 
!             if (volume_mesh%edge(ii,4) == -bb) then 
!                 write(11,'(I0,A,I0,A,I0)') 3,' ',volume_mesh%edge(ii,1) - 1,' ',volume_mesh%edge(ii,2) - 1
!             end if
!         end do 
!     end if 
! end do 




!Close file
close(11)

!Display
if (cm2dopt%dispt) then
    write(*,'(A)') '    {complete}'
end if
return 
end subroutine export_volume_mesh_SU2




!Import volume mesh subroutine ===========================
subroutine import_volume_mesh_flow(volume_mesh,cm2dopt)
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Local
integer(in) :: ii,vidx

!Open file 
open(11,file=cm2dopt%iopath//'grid')

!Read item quantities 
read(11,*) volume_mesh%ncell,volume_mesh%nedge,volume_mesh%nvtx

!Read mesh faces 
allocate(volume_mesh%edge(volume_mesh%nedge,4))
do ii=1,volume_mesh%nedge
    read(11,*) volume_mesh%edge(ii,:)
end do 

!Read mesh vertices 
allocate(volume_mesh%vertices(volume_mesh%nvtx,2))
do ii=1,volume_mesh%nvtx
    read(11,*) vidx,volume_mesh%vertices(ii,:)
end do 

!Read surface link data 
read(11,*) volume_mesh%nvtx_surf
allocate(volume_mesh%surf_vtx(volume_mesh%nvtx_surf))
allocate(volume_mesh%surf_vtx_seg(volume_mesh%nvtx_surf))
allocate(volume_mesh%surf_vtx_segfrac(volume_mesh%nvtx_surf))
do ii=1,volume_mesh%nvtx_surf 
    read(11,*) volume_mesh%surf_vtx(ii),volume_mesh%surf_vtx_seg(ii),volume_mesh%surf_vtx_segfrac(ii)
end do

!Close file 
close(11)

!Display
if (cm2dopt%dispt) then
    write(*,'(A)') '    {complete}'
end if
return 
end subroutine import_volume_mesh_flow




!Import flow gradients subroutine ===========================
subroutine import_flow_gradients(gradient_vol,volume_mesh,cm2dopt)
implicit none 

!Variables - Import
real(dp), dimension(:,:), allocatable :: gradient_vol
type(cm2d_options) :: cm2dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Local
integer(in) :: ii

!Allocate gradient array 
allocate(gradient_vol(volume_mesh%nvtx,2))

!Open file 
open(11,file=cm2dopt%iopath//'gradient.dat')

!Read gradients 
do ii=1,volume_mesh%nvtx
    read(11,*) gradient_vol(ii,:)
end do 

!Close file 
close(11)

!Display
if (cm2dopt%dispt) then
    write(*,'(A)') '    {complete}'
end if
return 
end subroutine import_flow_gradients




!Export surface gradients subroutine ===========================
subroutine export_surface_gradients(gradient_surf,surface_mesh,cm2dopt)
implicit none 

!Variables - Import
real(dp), dimension(:,:) :: gradient_surf
type(cm2d_options) :: cm2dopt
type(surface_data) :: surface_mesh

!Variables - Local
integer(in) :: ii

!Open file 
open(11,file=cm2dopt%iopath//'gradient_surf.dat')

!Write gradients 
do ii=1,surface_mesh%nvtx
    write(11,'(A,A,A)') real2F0_Xstring(gradient_surf(ii,1),12_in),' ',real2F0_Xstring(gradient_surf(ii,2),12_in)
end do 

!Close file 
close(11)

!Display
if (cm2dopt%dispt) then
    write(*,'(A)') '    {complete}'
end if
return 
end subroutine export_surface_gradients




!Write geometry check results subroutine ===========================
subroutine export_geometry_check(is_selfintersecting,cm2dopt)
implicit none 

!Variables - Import
logical :: is_selfintersecting
type(cm2d_options) :: cm2dopt

!Write
open(11,file=cm2dopt%iopath//'geometry_status') 
    if (is_selfintersecting) then 
        write(11,'(A,I0)') 'self intersection = ',1
    else
        write(11,'(A,I0)') 'self intersection = ',0
    end if 
close(11)
return 
end subroutine export_geometry_check




!Write status subroutine ===========================
subroutine export_status(cm2dopt)
implicit none 

!Variables - Import
type(cm2d_options) :: cm2dopt

!Write
open(11,file=cm2dopt%iopath//'cm2d_status') 
    write(11,'(I0)') cm2dopt%cm2dfailure
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
if (cm2dopt%dispt) then
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
if (cm2dopt%dispt) then
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

end module cellmesh2d_io_mod