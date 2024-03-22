!Quadtree Cutcell Surface Snapping Cell Containment Meshing Program (cell_mesh2d v2)
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 4.2
!Updated 22-03-2024

!Main
program cell_mesh2d
use cellmesh2d_mind_mod
use cellmesh2d_inflation_layer_mod
use cellmesh2d_mesh_generation_mod  
implicit none

!Variables
logical :: is_selfintersecting
type(cm2d_options) :: cm2dopt
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh
real(dp), dimension(:,:), allocatable :: gradient_vol,gradient_surf

!Set default paths 
call set_default_paths(cm2dopt)

!Set option defaults 
call set_default_options(cm2dopt)

!Get command arguments 
call get_process_arguments(cm2dopt)

!Import options
call cm2d_import_options(cm2dopt)

!Initialisation display
if (cm2dopt%dispt == 1) then
    write(*,'(A)') ' '
    write(*,'(A)')'+--------------------------------------------+'
    write(*,'(A)')'|              Cell Mesh 2D (v2)             |'
    write(*,'(A)')'|         2D Cut-Cell Mesh Generator         |'
    write(*,'(A)')'|        Version 0.9.6 || 22/03/2024         |'
    write(*,'(A)')'|                 Max Wood                   |'
    write(*,'(A)')'|           University of Bristol            |'
    write(*,'(A)')'|    Department of Aerospace Engineering     |'
    write(*,'(A)')'+--------------------------------------------+'
    write(*,'(A)') ' '
end if

!Process requested mode 
if (cm2dopt%mode == 'check') then !geometry check mode 

    !Display
    if (cm2dopt%dispt == 1) then
        write(*,'(A)') '         == geometry checking mode =='
    end if

    !Load object surface data
    if (cm2dopt%dispt == 1) then
        write(*,'(A)') '--> importing geometry data'
    end if
    call import_surface_geometry(surface_mesh,cm2dopt)

    !Preprocess surface geometry 
    if (cm2dopt%dispt == 1) then
        write(*,'(A)') '--> preprocessing geometry data'
    end if
    call preprocess_surface_mesh(surface_mesh,cm2dopt)

    !Check for self intersections 
    if (cm2dopt%dispt == 1) then
        write(*,'(A)') '--> testing for self intersections'
    end if
    is_selfintersecting = is_self_intersecting(surface_mesh,cm2dopt)
    if (cm2dopt%dispt == 1) then
        write(*,'(A,L,A)') '    {is self intersecting: ',is_selfintersecting,'}'
    end if
    
    !Export check results
    call export_geometry_check(is_selfintersecting,cm2dopt)
elseif (cm2dopt%mode == 'mesh') then !mesh generation mode

    !Display
    if (cm2dopt%dispt == 1) then
        write(*,'(A)') '         == mesh contruction mode =='
    end if

    !Import custom boundary condition zones
    if (cm2dopt%set_customBCs == 1) then 
        call cm2d_import_customBC_zones(cm2dopt)
    end if

    !Load object surface data
    if (cm2dopt%dispt == 1) then
        write(*,'(A)') '--> importing geometry data'
    end if
    call import_surface_geometry(surface_mesh,cm2dopt)

    !Preprocess surface mesh 
    if (cm2dopt%dispt == 1) then
        write(*,'(A)') '--> preprocessing surface mesh'
    end if
    call preprocess_surface_mesh(surface_mesh,cm2dopt)
    if (cm2dopt%cm2dfailure == 1) then
        call export_status(cm2dopt)
        stop 
    end if 

    !If inflation layer is requested then grow geometry to required height 
    if (cm2dopt%build_inflayer == 'yes') then 
        if (cm2dopt%dispt == 1) then
            write(*,'(A)') '--> inflating surface'
        end if
        call inflate_geometry(surface_mesh,cm2dopt)
    end if

    !Construct mesh 
    if (cm2dopt%dispt == 1) then
        write(*,'(A)') '--> constructing mesh'
    end if
    call cell_mesh2d_mesh(volume_mesh,surface_mesh,cm2dopt) !build base cutcell mesh 

    !Populate inflation layer mesh 
    if (cm2dopt%build_inflayer == 'yes') then 
        if (cm2dopt%dispt == 1) then
            write(*,'(A)') '--> meshing inflation layer'
        end if
        call mesh_inflation_layer(volume_mesh,surface_mesh,cm2dopt)
    end if 

    !Construct mesh cells if required 
    if ((cm2dopt%meshtype == 'su2_cutcell') .OR. (cm2dopt%meshtype == 'su2_dual')) then 
        call build_mesh_cells(volume_mesh)
    end if 

    !Build minD O mesh
    if (cm2dopt%meshtype == 'minD_O') then 
        call build_minD_Omesh(volume_mesh,surface_mesh,cm2dopt)
    end if 

    !Exit if critical failure detected 
    if (cm2dopt%cm2dfailure == 1) then 
        call export_status(cm2dopt)
        stop 
    end if 

    !Export items 
    call export_status(cm2dopt)
    call export_volume_mesh(volume_mesh,cm2dopt)
    if ((cm2dopt%meshtype == 'su2_cutcell') .OR. (cm2dopt%meshtype == 'su2_dual')) then 
        call export_volume_mesh_SU2(volume_mesh,cm2dopt)
    end if 
    if (cm2dopt%glink_con == 1) then 
        call export_vs2s_interpstruc(volume_mesh,surface_mesh,cm2dopt)
    end if 
elseif (cm2dopt%mode == 'project') then !volume to surface gradient projection mode 

    !Display
    if (cm2dopt%dispt == 1) then
        write(*,'(A)') '       == gradient projection mode =='
    end if

    !Load object surface data
    if (cm2dopt%dispt == 1) then
        write(*,'(A)') '--> importing geometry data'
    end if
    call import_surface_geometry(surface_mesh,cm2dopt)

    !Load volume mesh 
    if (cm2dopt%dispt == 1) then
        write(*,'(A)') '--> importing volume mesh'
    end if
    call import_volume_mesh_flow(volume_mesh,cm2dopt)

    !Load flow gradient file 
    if (cm2dopt%dispt == 1) then
        write(*,'(A)') '--> importing volume flow gradients'
    end if
    call import_flow_gradients(gradient_vol,volume_mesh,cm2dopt)

    !Preprocess surface mesh 
    if (cm2dopt%dispt == 1) then
        write(*,'(A)') '--> pre-processing the surface geometety'
    end if
    call orient_surface(surface_mesh,cm2dopt)
    call construct_surface_v2f(surface_mesh,cm2dopt)
    if (cm2dopt%cm2dfailure == 1) then
        call export_status(cm2dopt)
        stop 
    end if

    !Project volume gradients to the surface 
    if (cm2dopt%dispt == 1) then
        write(*,'(A)') '--> projecting volume flow gradients to the surface geometry'
    end if
    if (cm2dopt%glink_type == 'rbf') then 
        call project_gradients_RBF(gradient_surf,gradient_vol,volume_mesh,surface_mesh,4_in,&
        cm2dopt%ADTminNodedivsize,0.0d0,cm2dopt)
    elseif (cm2dopt%glink_type == 'int') then 
        call project_gradients_INT(gradient_surf,gradient_vol,volume_mesh,surface_mesh)
    end if 

    !Export surface projected gradients 
    if (cm2dopt%dispt == 1) then
        write(*,'(A)') '--> exporting surface gradients'
    end if
    call export_surface_gradients(gradient_surf,surface_mesh,cm2dopt)
else
    write(*,'(A,A)') '** unknown mode option requested : ',cm2dopt%mode
    stop
end if 

!Export status 
call export_status(cm2dopt)

!End
if (cm2dopt%dispt == 1) then
    stop
end if
end program cell_mesh2d