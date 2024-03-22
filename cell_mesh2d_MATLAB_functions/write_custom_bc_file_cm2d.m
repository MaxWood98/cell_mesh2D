%Function to write cell_mesh2d custom boundary conditions file 
%Max Wood - mw16116@bristol.ac.uk
%Univeristy of Bristol - Department of Aerospace Engineering

%Version 1.2
%Updated 22-03-2024

%Function -----------------------------------------------------------------
function [] = write_custom_bc_file_cm2d(BC_zones_loc,BC_zones_type)
    Nzone = size(BC_zones_type,1);
    fid = fopen(cm2dopt.bcondzone_filename,'w+');
    fprintf(fid,'%d \n',Nzone);
    for ii=1:Nzone
        fprintf(fid,'%d %f %f %f %f\n',BC_zones_type(ii),BC_zones_loc(ii,:));
    end
    fclose(fid);
end