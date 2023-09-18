!Quadtree Cutcell Surface Snapping Cell Containment Meshing Program (cell_mesh2d v2)
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 2.2
!Updated 18-09-2023

!Main
program cell_mesh
use cellmesh2d_mesh_generation_mod  
implicit none

!Variables
type(cm2d_options) :: cm2dopt
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh

!Get command arguments 
call get_process_arguments(cm2dopt)

!Import options
call cm2d_import_options(cm2dopt)

!Import custom boundary condition zones
if (cm2dopt%set_customBCs == 1) then 
    call cm2d_import_customBC_zones(cm2dopt)
end if

!Load object surface data
if (cm2dopt%dispt == 1) then
    write(*,*) '--> importing geometry data'
end if
call import_surface_geometry(surface_mesh,cm2dopt)

!Construct mesh
call cell_mesh2d_mesh(volume_mesh,surface_mesh,cm2dopt)

!Exit if failure detected 
if (cm2dopt%cm2dfailure == 1) then 
    call export_status(cm2dopt)
    stop 
end if 

!Export items 
call export_status(cm2dopt)
call export_volume_mesh(volume_mesh,cm2dopt)
if ((cm2dopt%meshfrmat == 'su2_cutcell') .OR. (cm2dopt%meshfrmat == 'su2_dual')) then 
    call export_volume_mesh_SU2(volume_mesh,cm2dopt)
end if 
if (cm2dopt%glink_con == 1) then 
    call export_vs2s_interpstruc(volume_mesh,surface_mesh,cm2dopt)
end if 

!End
if (cm2dopt%dispt == 1) then
    stop
end if
end program cell_mesh