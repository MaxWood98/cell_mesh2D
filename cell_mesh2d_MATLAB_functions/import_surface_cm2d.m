%Function to read cell_mesh2d surface file 
%Max Wood - mw16116@bristol.ac.uk
%Univeristy of Bristol - Department of Aerospace Engineering

%Version 1.2
%Updated 22-03-2024

%Function -----------------------------------------------------------------
function [Nvtx,Nedge,vertices,connectivity] = import_surface_cm2d(filename)

    %Read surface file
    fid = fopen(filename);
    surf = textscan(fid,'%f %f');
    fclose(fid);
    
    %Extract numver of vertices and objects
    Nvtx = surf{1}(1);
    Nedge = surf{2}(1);

    %Read vertices
    vertices = zeros(Nvtx,2);
    vertices(:,1) = surf{1}(2:Nvtx+1);
    vertices(:,2) = surf{2}(2:Nvtx+1);
    
    %Read connectivity
    connectivity = zeros(Nedge,2);
    connectivity(:,1) = surf{1}(Nvtx+2:Nvtx+Nedge+1);
    connectivity(:,2) = surf{2}(Nvtx+2:Nvtx+Nedge+1);
end