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
% cm2dop.meshfrmat = 'su2_dual';    %Mesh output format (flow2d / su2_cutcell / su2_dual)
cm2dop.meshfrmat = 'flow2d';
% cm2dop.meshfrmat = 'su2_cutcell';

%Cut-Cell mesh options ====================================================
%Quadtree options
cm2dop.nrefine = 12;              %Maximum refinement level 
cm2dop.nrefineB = 2;              %Maximum additional refinement levels in high curvature regions
cm2dop.ncell_max = 200000;        %Maximum number of cells
cm2dop.nrflood_i = 6;             %Refinement adjacency flooding iterations at the first refinement
cm2dop.nrflood_f = 4;             %Refinement adjacency flooding iterations at the final refinement
cm2dop.nrflood_b = 3;             %Refinement adjacency flooding iterations on boosted refinement
cm2dop.fbound = 15;               %Far field distance from object centre  
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
cm2dop.glink_con = 1;             %Construct volume to surface gradient interpolation (1 = yes | 0 = no)
cm2dop.glink_nnn = 10;            %Number of nearest neighbours to use for volume to surface gradient interpolation
cm2dop.glink_nsmooth = 4;         %Number of vertices each side used to smooth the gradient at each surface vertex

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
system('cell_mesh2d');

%% Load mesh

%Load mesh
filepath = 'io/grid';
[Ncell,Nedge,Nvtx,edge,vtx,cell_lr] = import_mesh_cm2d(filepath);

%Load gradient link structure
% [vs2s_interp] = import_vs2s_interp_cm2d('io/vs2s_interp');

% edge(edge == 0) = nan;
% for ii=1:Nedge
%     if cell_lr(ii,1) == -1 
%         edge(ii,:) = nan;
%     end
% end 

% %Check cells
% for ii=1:Nedge
%     if cell_lr(ii,1) == 0 || cell_lr(ii,2) == 0
%         disp('** zero adjacent edge')
%         cell_lr(ii,:)
%     end
% end 

%% Plot mesh

%Setup figure
cla reset 
hold on

%Plot mesh
patch('vertices',vtx,'faces',edge,'EdgeAlpha',1.0,'Marker','none');
% plot(vtx(:,1),vtx(:,2),'r.')

%Read input surface file 
[~,~,vertices,connectivity] = import_cell_mesh2d_surface('io/cell_mesh2d_surface.dat');

%Plot object surface 
patch('vertices',vertices,'faces',connectivity,'EdgeAlpha',0.5,'Marker','none','EdgeColor',[0.1 0.1 1],'MarkerEdgeColor','b');
% plot(vertices(:,1),vertices(:,2),'b.')
% plot(vertices(360,1),vertices(360,2),'b.','markersize',20)

%Plot boundary conditions 
for ii=1:Nedge
    if cell_lr(ii,1) == -1 %wall
        patch('vertices',vtx,'faces',edge(ii,:),'EdgeAlpha',1.0,'Marker','none','Edgecolor','c');
    elseif cell_lr(ii,1) == -2 %far field
        patch('vertices',vtx,'faces',edge(ii,:),'EdgeAlpha',1.0,'Marker','none','Edgecolor','b');
    elseif cell_lr(ii,1) == -3 || cell_lr(ii,1) == -5 || cell_lr(ii,1) == -6 %inflow
        patch('vertices',vtx,'faces',edge(ii,:),'EdgeAlpha',1.0,'Marker','none','Edgecolor','g');
    elseif cell_lr(ii,1) == -4 || cell_lr(ii,1) == -7 %outflow
        patch('vertices',vtx,'faces',edge(ii,:),'EdgeAlpha',1.0,'Marker','none','Edgecolor','r');
    end
end

%Debug plots ==============================================================


% vrcurv = load('io/rcurv.dat');

% celloftype = load('io/celloftype');
% plot(celloftype(:,1),celloftype(:,2),'r.')

% vtxexternal = load('io/vtxexternal');
% plot(vtxexternal(:,1),vtxexternal(:,2),'r.')

% surfvtx = load('io/surfvtx');
% plot(surfvtx(:,1),surfvtx(:,2),'r.')



% cmid = zeros(Ncell,2);
% cmidls = zeros(Ncell,1);
% for ii=1:Nedge
%     cl = cell_lr(ii,1);
%     cr = cell_lr(ii,2);
%     v1 = edge(ii,1);
%     v2 = edge(ii,2);
%     emidx = 0.5*(vtx(v1,1) + vtx(v2,1));
%     emidy = 0.5*(vtx(v1,2) + vtx(v2,2));
%     dx = vtx(v2,1) - vtx(v1,1);
%     dy = vtx(v2,2) - vtx(v1,2);
%     ledge = sqrt(dx^2 + dy^2);
%     if cl > 0
%         cmid(cl,1) = cmid(cl,1) + emidx*ledge;
%         cmid(cl,2) = cmid(cl,2) + emidy*ledge;
%         cmidls(cl) = cmidls(cl) + ledge;
%     end
%     if cr > 0
%         cmid(cr,1) = cmid(cr,1) + emidx*ledge;
%         cmid(cr,2) = cmid(cr,2) + emidy*ledge;
%         cmidls(cr) = cmidls(cr) + ledge;
%     end
% end
% cmid(:,:) = cmid(:,:)./cmidls(:);
% for ii=1:Nedge
%     etgt = ii;
%     if cell_lr(etgt,1) == 6479 || cell_lr(etgt,2) == 6479
%         v1 = edge(etgt,1);
%         v2 = edge(etgt,2);
%         emidx = 0.5*(vtx(v1,1) + vtx(v2,1));
%         emidy = 0.5*(vtx(v1,2) + vtx(v2,2));
%         cadj = cell_lr(etgt,2);
% 
%         ledge = norm(vtx(v1,:) - vtx(v2,:));
% 
%         if cadj > 0
%             plot([emidx cmid(cadj,1)],[emidy cmid(cadj,2)],'b','linewidth',2)
%         else
%             plot([vtx(v1,1) vtx(v2,1)],[vtx(v1,2) vtx(v2,2)],'r','linewidth',2)
%         end
%     end
% end


% cell_remove = load('io/cell_remove');
% for ii=1:Nedge
%     etgt = ii;
%     if cell_remove(cell_lr(etgt,2)) == 1 && cell_lr(etgt,1) < 0 
%         v1 = edge(etgt,1);
%         v2 = edge(etgt,2);
%         plot([vtx(v1,1) vtx(v2,1)],[vtx(v1,2) vtx(v2,2)],'r','linewidth',2)
%     end
% end

% ctgt = 6745;
% for ii=1:Nedge
%     etgt = ii;
%     if cell_lr(etgt,1) == ctgt || cell_lr(etgt,2) == ctgt 
%         v1 = edge(etgt,1);
%         v2 = edge(etgt,2);
%         plot([vtx(v1,1) vtx(v2,1)],[vtx(v1,2) vtx(v2,2)],'r','linewidth',2)
%     end
% end

% ctgt = 12;
% for ii=1:Nedge
%     etgt = ii;
%     if cell_lr(etgt,1) == ctgt || cell_lr(etgt,2) == ctgt 
%         v1 = edge(etgt,1);
%         v2 = edge(etgt,2);
%         plot([vtx(v1,1) vtx(v2,1)],[vtx(v1,2) vtx(v2,2)],'r','linewidth',2)
%     end
% end


% vtgt = 17493;
% plot(vtx(vtgt,1),vtx(vtgt,2),'r.','markersize',20)
% 
% vtgt = 17397;
% plot(vtx(vtgt,1),vtx(vtgt,2),'r.','markersize',20)
% 
% vtgt = 17365;
% plot(vtx(vtgt,1),vtx(vtgt,2),'r.','markersize',20)
% 
% vtgt = 17889;
% plot(vtx(vtgt,1),vtx(vtgt,2),'r.','markersize',20)
% 
% vtgt = 17898;
% plot(vtx(vtgt,1),vtx(vtgt,2),'r.','markersize',20)
% 
% vtgt = 18198;
% plot(vtx(vtgt,1),vtx(vtgt,2),'r.','markersize',20)
% 
% vtgt = 18303;
% plot(vtx(vtgt,1),vtx(vtgt,2),'r.','markersize',20)
% 
% vtgt = 18816;
% plot(vtx(vtgt,1),vtx(vtgt,2),'r.','markersize',20)


% plot(0.99999999999511724   ,     1.3427734375000000E-003 ,'bo')
% plot(0.99999999999511724   ,     1.2207031250000000E-003,'bo')
% patch('vertices',vertices,'faces',connectivity(131,:),'EdgeAlpha',1,'Marker','none','EdgeColor',[0.1 0.1 1],'MarkerEdgeColor','b');


% % 3086  0.89072328226249309       -1.5625000000000000E-002
% plot(0.89072328226249309    ,   -1.5625000000000000E-002,'r*')
% 
% plot(0.90625000000000000    ,   -1.5625000000000000E-002,'bo')
% plot(0.89062500000000000    ,   -1.5625000000000000E-002,'go')
% 
% plot(cmid(1466,1),cmid(1466,2),'g*')
% plot(cmid(1465,1),cmid(1465,2),'m*')

% [vs2s_interp] = import_vs2s_interp_cm2d('io/vs2s_interp');
% vtgt = 269;
% volpoints = vs2s_interp.idata{vtgt,1}.volpoints(:);
% plot(vtx(volpoints,1),vtx(volpoints,2),'r.')
% plot(vertices(vtgt,1),vertices(vtgt,2),'b*')



% interp_curve = load('io/interp_curve.dat');
% plot(interp_curve(:,1),interp_curve(:,2),'r','linewidth',2)
% 
% vbase = load('io/interp_curve_vbase.dat');
% plot(vertices(vbase(:),1),vertices(vbase(:),2),'g.','markersize',15)

% vtgt = 21;
% plot(vertices(vbase(vtgt),1),vertices(vbase(vtgt),2),'b.','markersize',15)

% ctgt = 3062;
% for ii=1:Nedge
%     etgt = ii;
%     if cell_lr(etgt,1) == ctgt || cell_lr(etgt,2) == ctgt 
%         etgt
%         v1 = edge(etgt,1);
%         v2 = edge(etgt,2);
%         plot([vtx(v1,1) vtx(v2,1)],[vtx(v1,2) vtx(v2,2)],'r','linewidth',2)
%     end
% end


% verticesT = load('io/vtxtest.dat');
% plot(verticesT(:,1),verticesT(:,2),'r.')
% 



% for ii=1:Nedge
%     etgt = ii;
%     if cell_lr(etgt,1) < 0 
%         v1 = edge(etgt,1);
%         v2 = edge(etgt,2);
%         emidx = 0.5*(vtx(v1,1) + vtx(v2,1));
%         emidy = 0.5*(vtx(v1,2) + vtx(v2,2));
%         dx = vtx(v2,1) - vtx(v1,1);
%         dy = vtx(v2,2) - vtx(v1,2);
%         plot([emidx emidx+dy],[emidy emidy-dx],'r','linewidth',2)
%     end
% end


%Debug plots ==============================================================


%Format
axis equal
axis tight
xlabel('x')
ylabel('y')
hold off

% axis([-0.9789   -0.9747   -0.0020    0.0022]);


