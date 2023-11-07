%Quadtree Cutcell Meshing Program (Cell Mesh 2D) Interface
%Max Wood - mw16116@bristol.ac.uk
%Univeristy of Bristol - Department of Aerospace Engineering

%Version 8.0
%Updated 21-07-2023

%Reset workspace
clearvars
clc

%TODO
%- mind mesh add surface zone construction on far field regions 

%% Input

%General options 
cm2dop.condisp = 1;               %Toggle console display (1 = yes || 0 = no)

%Mesh format options 
cm2dop.meshtype = 0;              %Type of mesh (0 = cutcell | 1 = minD block mesh)
cm2dop.meshinout = 'out';         %Mesh inside or outside of geometry (default out)
cm2dop.surface_dir = 'in';        %Surface normal direction switch in to / out of the mesh domain (default in)
cm2dop.boundary_dir = 'in';       %Boundary normal direction switch in to / out of the mesh domain (default in)
cm2dop.meshfrmat = 'cutcell';    %Mesh output format (cutcell / su2_cutcell / su2_dual)
% cm2dop.meshfrmat = 'su2_dual';
% cm2dop.meshfrmat = 'su2_cutcell';

%Cut-Cell mesh options ====================================================
%Quadtree options
cm2dop.nrefine = 12;              %Maximum refinement level 
cm2dop.nrefineB = 2;              %Maximum additional refinement levels in high curvature regions
cm2dop.ncell_max = 200000;        %Maximum number of cells
cm2dop.nrflood_i = 6;             %Refinement adjacency flooding iterations at the first refinement
cm2dop.nrflood_f = 4;             %Refinement adjacency flooding iterations at the final refinement
cm2dop.nrflood_b = 2;             %Refinement adjacency flooding iterations on boosted refinement
cm2dop.fbound = 12;               %Far field distance from object centre  
cm2dop.coffset = [0.0 0.0];       %Object/mesh centre offset (x / y)

%Custom domain bound specification 
cm2dop.set_mbounds = 0;           %Toggle forcing of custom mesh domain bounds
cm2dop.xmin = 0.2;                %xmin
cm2dop.xmax = 11.8;               %xmax
cm2dop.ymin = 0.2;                %ymin
cm2dop.ymax = 11.8;               %ymin

%Mesh cleaning options
cm2dop.eminlen = 1e-8;            %Minimum edge length as a fraction of an undeformed cell edge length at each refienemnt level 
cm2dop.cminvol = 0.01;            %Volume fraction of an undeformed cell at each refinement level below which a cell is classed as a sliver cell

%Mesh geometry intersection options
cm2dop.enintmax = 50;             %Maximum number of mesh-geometry intersections on each mesh edge
cm2dop.eintpad = 0.0;             %Edge-geometry intersection search zone padding as a fraction of the edge length 
cm2dop.int_coin_tol = 1e-8;       %Intersection co-incidence tollerance as fraction of item length  

%Surface format options
cm2dop.surftype = 0;              %Geometry surface type (0 = 'simplified' | 1 = 'exact') 
cm2dop.surfRcurvM = 1.0;          %Surface curvature multiplier
cm2dop.surfRcurvNpts = 10;        %Number of vertices used to estimate local surface curvature

%Mesh smoothing options
cm2dop.nsstype = 0;               %Near surface smoothing type (0 = 'none' | 1 = 'Laplacian')
cm2dop.nsvtxflood = 2;            %Vertex selection flooding iterations from surfaces
cm2dop.nlpsmooth = 3;             %Smoothing iterations 

%ADtree options
cm2dop.adtree_spad = 0.0;         %Maximum padding size of adtree search bounding boxes as multiple of cell edge length
cm2dop.adtree_maxd = 4;           %AD tree maximum depth in terms of dimension cycles (tree is 4d)

%Gradient linking options
cm2dop.glink_con = 0;             %Construct volume to surface gradient interpolation (1 = yes | 0 = no)
cm2dop.glink_nnn = 10;            %Number of nearest neighbours to use for volume to surface gradient interpolation
cm2dop.glink_nsmooth = 0;         %Number of vertices each side used to smooth the gradient at each surface vertex
cm2dop.glink_RBF_relax = 1e-6;    %RBF interpolation smoothing relaxation parameter

%Boundary condition options ===============================================
%Boundary condition options 
cm2dop.set_custom_bc = 0;         %Set custom boundary conditions in specifed regions (1 = yes | 0 = no)
cm2dop.rem_ffzones = 0;           %Remove any region of the mesh connected to a far field boundary condition
cm2dop.rem_nczones = 0;           %Remove any region of the mesh connected to a non custom set boundary condition
cm2dop.rem_iszones = 0;           %Remove any isolated region of the mesh connected only to a wall boundary condition 

%Boundary condition zone bounds [xmin xmax ymin ymax]
BC_zones_loc = [0.0 1.0 0.5 3.5;
                5.5 8.5 11.5 12.5];

%Boundary condition zone conditions
% -1 = wall
% -2 = far field
% -3 = mass flux inflow 
% -4 = mass flux outflow
% -5 = stagnation state inflow
% -6 = freestream inflow
% -7 = back pressure outflow 
BC_zones_type = [-2 ; -6];



%% Meshing function

%Write options files
write_input_file_cm2d(cm2dop);
if cm2dop.set_custom_bc == 1
    write_custom_bc_file_cm2d(BC_zones_loc,BC_zones_type);
end 

%Call meshing function 
system('cell_mesh2d mesh');
system('cell_mesh2d project');

%% Load mesh

%Load mesh
filepath = 'io/grid';
[Ncell,Nedge,Nvtx,edge,vtx,cell_lr] = import_mesh_cm2d(filepath);

%Load gradient link structure
% [vs2s_interp] = import_vs2s_interp_cm2d('io/vs2s_interp');

%% Plot mesh

%Setup figure
cla reset 
hold on

%Plot mesh
patch('vertices',vtx,'faces',edge,'EdgeAlpha',1.0,'Marker','none');

%Read input surface file 
[~,~,vertices,connectivity] = import_cell_mesh2d_surface('io/cell_mesh2d_surface.dat');

%Plot object surface 
% patch('vertices',vertices,'faces',connectivity,'EdgeAlpha',0.5,'Marker','none','EdgeColor',[0.1 0.1 1],'MarkerEdgeColor','b');

%Plot boundary conditions 
for ii=1:Nedge
    if cell_lr(ii,1) == -1 %wall
        % patch('vertices',vtx,'faces',edge(ii,:),'EdgeAlpha',1.0,'Marker','none','Edgecolor','c');
    elseif cell_lr(ii,1) == -2 %far field
        patch('vertices',vtx,'faces',edge(ii,:),'EdgeAlpha',1.0,'Marker','none','Edgecolor','b');
    elseif cell_lr(ii,1) == -3 || cell_lr(ii,1) == -5 || cell_lr(ii,1) == -6 %inflow
        patch('vertices',vtx,'faces',edge(ii,:),'EdgeAlpha',1.0,'Marker','none','Edgecolor','g');
    elseif cell_lr(ii,1) == -4 || cell_lr(ii,1) == -7 %outflow
        patch('vertices',vtx,'faces',edge(ii,:),'EdgeAlpha',1.0,'Marker','none','Edgecolor','r');
    end
end

% %Plot surface and volume gradients 
grad_surf = load('io/gradient_surf.dat');
quiver(vertices(:,1),vertices(:,2),grad_surf(:,1),grad_surf(:,2),0,'r','maxheadsize',0.01)
grad_vol = load('io/gradient.dat');
quiver(vtx(:,1),vtx(:,2),grad_vol(:,1),grad_vol(:,2),0,'b','maxheadsize',0.01)


%Format
axis equal
axis tight
xlabel('x')
ylabel('y')
hold off

axis([-0.1699    1.1425   -0.6533    0.6591]);