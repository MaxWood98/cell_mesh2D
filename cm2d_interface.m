%Quadtree Cutcell Meshing Program (Cell Mesh 2D) Interface
%Max Wood - mw16116@bristol.ac.uk
%Univeristy of Bristol - Department of Aerospace Engineering

%Version 8.0
%Updated 21-07-2023

%Reset workspace
clearvars
clc

%TODO
%- add input for cell aspect ratio merge denial
%- mind mesh add surface zone construction on far field regions 

%% Input

%General options 
cm2dopt.condisp = 1;               %Toggle console display (1 = yes || 0 = no)

%Mesh format options 
cm2dopt.meshtype = 0;              %Type of mesh (0 = cutcell | 1 = minD block mesh)
cm2dopt.meshinout = 'out';         %Mesh inside or outside of geometry (default out)
cm2dopt.surface_dir = 'in';        %Surface normal direction switch in to / out of the mesh domain (default in)
cm2dopt.boundary_dir = 'in';       %Boundary normal direction switch in to / out of the mesh domain (default in)
cm2dopt.meshfrmat = 'cutcell';     %Mesh output format (cutcell / su2_cutcell / su2_dual)
% cm2dopt.meshfrmat = 'su2_dual';
% cm2dopt.meshfrmat = 'su2_cutcell';

%Cut-Cell mesh options ====================================================
%Quadtree options
cm2dopt.nrefine = 11;              %Maximum refinement level 
cm2dopt.nrefineB = 0;              %Maximum additional refinement levels in high curvature regions
cm2dopt.ncell_max = 200000;        %Maximum number of cells
cm2dopt.nrflood_i = 8;             %Refinement adjacency flooding iterations at the first refinement
cm2dopt.nrflood_f = 10;            %Refinement adjacency flooding iterations at the final refinement
cm2dopt.nrflood_b = 1;             %Refinement adjacency flooding iterations on boosted refinement
cm2dopt.fbound = 15;               %Far field distance from object centre  
cm2dopt.coffset = [0.0 0.0];       %Object/mesh centre offset (x / y)

%Custom domain bound specification 
cm2dopt.set_mbounds = 0;           %Toggle forcing of custom mesh domain bounds
cm2dopt.xmin = 0.2;                %xmin
cm2dopt.xmax = 11.8;               %xmax
cm2dopt.ymin = 0.2;                %ymin
cm2dopt.ymax = 11.8;               %ymin

%Mesh cleaning options
cm2dopt.eminlen = 1e-8;            %Minimum edge length as a fraction of an undeformed cell edge length at each refienemnt level 
cm2dopt.cminvol = 0.1;             %Volume fraction of an undeformed cell at each refinement level below which a cell is classed as a sliver cell

%Mesh geometry intersection options
cm2dopt.enintmax = 50;             %Maximum number of mesh-geometry intersections on each mesh edge
cm2dopt.eintpad = 0.0;             %Edge-geometry intersection search zone padding as a fraction of the edge length 
cm2dopt.int_coin_tol = 1e-8;       %Intersection co-incidence tollerance as fraction of item length  

%Surface format options
cm2dopt.surftype = 0;              %Geometry surface type (0 = 'simplified' | 1 = 'exact') 
cm2dopt.surfRcurvM = 1.0;          %Surface curvature multiplier
cm2dopt.surfRcurvNpts = 10;        %Number of vertices used to estimate local surface curvature

%Mesh smoothing options
cm2dopt.nsstype = 0;               %Near surface smoothing type (0 = 'none' | 1 = 'Laplacian')
cm2dopt.nsvtxflood = 2;            %Vertex selection flooding iterations from surfaces
cm2dopt.nlpsmooth = 3;             %Smoothing iterations 

%ADtree options
cm2dopt.adtree_spad = 0.0;         %Maximum padding size of adtree search bounding boxes as multiple of cell edge length
cm2dopt.adtree_maxd = 10;          %AD tree maximum depth in terms of dimension cycles (tree is 4d)
cm2dopt.adtree_mindivnsize = 10;   %AD tree minimum divisible node size 

%Gradient linking options
cm2dopt.glink_con = 0;             %Construct volume to surface gradient interpolation (1 = yes | 0 = no)
cm2dopt.glinktype = 'rbf';         %Gradient linking type (rbf or int)
cm2dopt.glink_nnn = 10;            %Number of nearest neighbours to use for RBF volume to surface gradient interpolation
cm2dopt.glink_nsmooth = 0;         %Number of vertices each side used to smooth the gradient at each surface vertex

%Radial Basis function Options 
cm2dopt.RBF_rsup = 50.0;           %RBF interpolation support radius as a multiple of the maximum distance between points
cm2dopt.RBF_relaxD = 0.05;         %RBF interpolation relaxation distance as fraction of the maximum distance between points
cm2dopt.RBF_relaxP = 0.5;          %RBF interpolation relaxation parameter

%Boundary condition options ===============================================
%Boundary condition options 
cm2dopt.set_custom_bc = 0;         %Set custom boundary conditions in specifed regions (1 = yes | 0 = no)
cm2dopt.rem_ffzones = 0;           %Remove any region of the mesh connected to a far field boundary condition
cm2dopt.rem_nczones = 0;           %Remove any region of the mesh connected to a non custom set boundary condition
cm2dopt.rem_iszones = 0;           %Remove any isolated region of the mesh connected only to a wall boundary condition 

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
write_input_file_cm2d(cm2dopt);
if cm2dopt.set_custom_bc == 1
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
patch('vertices',vertices,'faces',connectivity,'EdgeAlpha',0.5,'Marker','.','EdgeColor',[0.1 0.1 1],'MarkerEdgeColor','b');

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

%Plot surface and volume gradients 
grad_surf = load('io/gradient_surf.dat');
quiver(vertices(:,1),vertices(:,2),grad_surf(:,1),grad_surf(:,2),0,'r','maxheadsize',0.01)
grad_vol = load('io/gradient.dat');
quiver(vtx(:,1),vtx(:,2),grad_vol(:,1),grad_vol(:,2),0,'b','maxheadsize',0.01)

% %Plot cell
% % ctgt = 4475;
% % ctgt = 4654;
% ctgt = 4653;
% for ii=1:Nedge
%     if cell_lr(ii,1) == ctgt || cell_lr(ii,2) == ctgt
%         % edge(ii,:)
%         % cell_lr(ii,:)
%         patch('vertices',vtx,'faces',edge(ii,:),'EdgeAlpha',1.0,'Marker','.','Edgecolor','g');
%     end
% end 


%Format
axis equal
axis tight
xlabel('x')
ylabel('y')
hold off

axis([-0.1699    1.1425   -0.6533    0.6591]);