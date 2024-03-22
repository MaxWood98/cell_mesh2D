%cell_mesh2d run function
%Max Wood - mw16116@bristol.ac.uk
%Univeristy of Bristol - Department of Aerospace Engineering

%Version 0.1
%Updated 22-03-2024

%Function -----------------------------------------------------------------
function [] = cell_mesh2d_run(mode,options_filepath,surface_filepath)
    
    %Construct command
    command = ['cell_mesh2d ',mode,' -o ',options_filepath,' -s ',surface_filepath];

    %Mesh   
    system(command);
end