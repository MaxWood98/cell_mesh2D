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

%Filepaths
cm2dopt.options_filepath = 'cell_mesh2d_options';                 %Options filepath
cm2dopt.surface_filepath = 'io/cell_mesh2d_surface_SHOPT.dat';     %Surface filepath
cm2dopt.bcondzone_filename = 'cell_mesh2d_bcond_zones';           %Boundary conditions filepath

%General options 
cm2dopt.condisp = 'yes';           %Toggle console display (yes | no)

%Mesh format options 
cm2dopt.meshtype = 'cutcell';      %Type of mesh (cutcell | su2_cutcell | su2_dual | minD_O)
cm2dopt.meshinout = 'out';         %Mesh inside or outside of geometry (default out)
cm2dopt.surface_dir = 'in';        %Surface normal direction switch in to / out of the mesh domain (default in)
cm2dopt.boundary_dir = 'in';       %Boundary normal direction switch in to / out of the mesh domain (default in)

%Cut-Cell mesh options ====================================================
%Quadtree options
cm2dopt.nrefine = 12;              %Maximum refinement level 
cm2dopt.nrefineB = 2;              %Maximum additional refinement levels in high curvature regions
cm2dopt.ncell_max = 200000;        %Maximum number of cells
cm2dopt.nrflood_i = 7;             %Refinement adjacency flooding iterations at the first refinement
cm2dopt.nrflood_f = 7;             %Refinement adjacency flooding iterations at the final refinement
cm2dopt.nrflood_b = 4;             %Refinement adjacency flooding iterations on boosted refinement
cm2dopt.fbound = 20;               %Far field distance from object centre  
cm2dopt.coffset = [0.0 0.0];       %Object/mesh centre offset (x y)

%Custom domain bound specification 
cm2dopt.set_mbounds = 'no';        %Set custom mesh domain bounds (yes | no)
cm2dopt.xmin = 0.2;                %xmin
cm2dopt.xmax = 11.8;               %xmax
cm2dopt.ymin = 0.2;                %ymin
cm2dopt.ymax = 11.8;               %ymin

%Mesh cleaning options
cm2dopt.eminlen = 1e-8;            %Minimum edge length as a fraction of an undeformed cell edge length at each refienemnt level 
cm2dopt.srfeminlen = 0.005;        %Minimum surface edge length as a fraction of an undeformed cell edge length at each refienemnt level 
cm2dopt.cminvol = 0.1;             %Volume fraction of an undeformed cell at each refinement level below which a cell is classed as a sliver cell
cm2dopt.srfinclean = 'no';         %Allow cleaning of the input surface geometry 

%Mesh geometry intersection options
cm2dopt.enintmax = 50;             %Maximum number of mesh-geometry intersections on each mesh edge
cm2dopt.eintpad = 0.0;             %Edge-geometry intersection search zone padding as a fraction of the edge length 
cm2dopt.int_coin_tol = 1e-8;       %Intersection co-incidence tollerance as fraction of item length  

%Surface format options
cm2dopt.surftype = 'simplified';   %Geometry surface type (simplified | exact) 
cm2dopt.surfRcurvM = 1.0;          %Surface curvature multiplier
cm2dopt.surfRcurvNpts = 10;        %Number of vertices used to estimate local surface curvature

%Inflation layer options 
cm2dopt.build_inflayer = 'no';     %Toggle construction of an inflation layer (yes | no)
cm2dopt.inflayer_height = 0.04;    %Inflation layer approximate height
cm2dopt.inflayer_nlayer = 0;       %Number of layers within the inflation layer (determine automatically if zero)
cm2dopt.inflayer_wd = 0.1;         %Inflation layer differencing weight (0->1)
cm2dopt.inflayer_cvxdp = 0.1;      %Inflation layer convex differencing penalty (0->1)
cm2dopt.inflayer_we = 0.5;         %Inflation layer evening weight (0->1)
cm2dopt.inflayer_ebcbase = 0.1;    %Inflation layer concave evening bias base value
cm2dopt.inflayer_enflood = 10;     %Inflation layer evening number of flood iterations
cm2dopt.inflayer_ensubiter = 10;   %Inflation layer evening number of sub-iterations
cm2dopt.inflayer_lreb = 1.0;       %Inflation layer evening length ratio bound (>=1)

%Mesh smoothing options
cm2dopt.nsstype = 'none';          %Near surface smoothing type (none | laplacian)
cm2dopt.nsvtxflood = 5;            %Vertex selection flooding iterations from surfaces
cm2dopt.nlpsmooth = 10;            %Smoothing iterations 

%ADtree options
cm2dopt.adtree_spad = 0.0;         %Maximum padding size of adtree search bounding boxes as multiple of cell edge length
cm2dopt.adtree_maxd = 10;          %AD tree maximum depth in terms of dimension cycles (tree is 4d)
cm2dopt.adtree_mindivnsize = 10;   %AD tree minimum divisible node size 

%Gradient projection options
cm2dopt.glink_con = 'no';          %Construct and export volume to surface gradient interpolation (yes | no)
cm2dopt.glinktype = 'int';         %Gradient linking type (rbf | int)
cm2dopt.glink_nnn = 80;            %Number of nearest neighbours to use for RBF volume to surface gradient interpolation
cm2dopt.glink_nsmooth = 0;         %Number of vertices each side used to smooth the gradient at each surface vertex

%Radial Basis function Options 
cm2dopt.RBF_rsup = 50.0;           %RBF interpolation support radius as a multiple of the maximum distance between points
cm2dopt.RBF_relaxD = 0.01;         %RBF interpolation relaxation distance as fraction of the maximum distance between points
cm2dopt.RBF_relaxP = 0.01;         %RBF interpolation relaxation parameter

%Boundary condition options ===============================================
%Boundary condition options 
cm2dopt.set_custom_bc = 'no';      %Set custom boundary conditions in specifed regions (yes | no)
cm2dopt.rem_ffzones = 'no';        %Remove any region of the mesh connected to a far field boundary condition (yes | no)
cm2dopt.rem_nczones = 'no';        %Remove any region of the mesh connected to a non custom set boundary condition (yes | no)
cm2dopt.rem_iszones = 'no';        %Remove any isolated region of the mesh connected only to a wall boundary condition (yes | no)

%Boundary condition zone bounds [xmin xmax ymin ymax]
% BC_zones_loc = [-1.0 1.0 0.5 3.5;
%                 5.5 8.5 11.5 12.5];

% BC_zones_loc = [-1.0 1.0 0.5 7.5;
%                 5.5 8.5 -1.0 9.0];

%Boundary condition zone conditions
% -1 = wall
% -2 = far field
% -3 = mass flux inflow 
% -4 = mass flux outflow
% -5 = stagnation state inflow
% -6 = freestream inflow
% -7 = back pressure outflow 
BC_zones_type = [-6 ; -7];



%% Meshing function

%Write options files
write_input_file_cm2d(cm2dopt);
if cm2dopt.set_custom_bc == 1
    write_custom_bc_file_cm2d(BC_zones_loc,BC_zones_type);
end 

%Call meshing function 
cell_mesh2d_run('mesh',cm2dopt.options_filepath,cm2dopt.surface_filepath);
% cell_mesh2d_run('project',cm2dopt.options_filepath,cm2dopt.surface_filepath);

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

% vtxd = load('io/vtxd');

%Plot mesh
patch('vertices',vtx,'faces',edge,'EdgeAlpha',1.0,'Marker','none');
% patch('vertices',[vtx vtxd],'faces',edge,'EdgeAlpha',1.0,'Marker','none');



% vtgt = 366+1;
% plot(vtx(vtgt,1),vtx(vtgt,2),'r*');
% 
% vtgt = 368+1;
% plot(vtx(vtgt,1),vtx(vtgt,2),'b*');
% 
% vtgt = 788+1;
% plot(vtx(vtgt,1),vtx(vtgt,2),'r*');
% 
% vtgt = 797+1;
% plot(vtx(vtgt,1),vtx(vtgt,2),'b*');
% 
% vtgt = 487+1;
% plot(vtx(vtgt,1),vtx(vtgt,2),'r*');
% 
% vtgt = 367+1;
% plot(vtx(vtgt,1),vtx(vtgt,2),'b*');


% vtgt = 5115;
% plot(vtx(vtgt,1),vtx(vtgt,2),'b*');
% 
% vtgt = 5116;
% plot(vtx(vtgt,1),vtx(vtgt,2),'g*');
% 
% vtgt = 1;
% plot(vtx(vtgt,1),vtx(vtgt,2),'r*');


% vtgt = 5115;
% plot(vtx(vtgt,1),vtx(vtgt,2),'b*');
% 
% vtgt = 1;
% plot(vtx(vtgt,1),vtx(vtgt,2),'g*');
% 
% vtgt = 2;
% plot(vtx(vtgt,1),vtx(vtgt,2),'r*');
% 
% vtgt = 5117;
% plot(vtx(vtgt,1),vtx(vtgt,2),'k*');


%Read input surface file 
[~,~,vertices,connectivity] = import_surface_cm2d(cm2dopt.surface_filepath);

%Plot object surface 
patch('vertices',vertices,'faces',connectivity,'EdgeAlpha',0.5,'Marker','.','EdgeColor',[0.1 0.1 1],'MarkerEdgeColor','b');

% vtgt = 1;
% plot(vertices(vtgt,1),vertices(vtgt,2),'g*');
% 
% vtgt = 3;
% plot(vertices(vtgt,1),vertices(vtgt,2),'g*');


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
% % grad_surf = load('io/gradient_surf.dat');
% % quiver(vertices(:,1),vertices(:,2),grad_surf(:,1),grad_surf(:,2),0,'r','maxheadsize',0.01)
% grad_vol = load('io/gradient.dat');
% grad_valid = zeros(length(grad_vol(:,1)),1);
% for ii=1:Nedge
%     if cell_lr(ii,1) == -1 %wall
%         grad_valid(edge(ii,:)) = 1;
%     end
% end
% for ii=1:length(grad_vol(:,1))
%     if grad_valid(ii) == 0
%         grad_vol(ii,:) = 0;
%     end
% end
% quiver(vtx(:,1),vtx(:,2),grad_vol(:,1),grad_vol(:,2),0.5,'b','maxheadsize',0.01)




% %Plot cell
% % ctgt = 4475;
% % ctgt = 4654;
% ctgt = 14770;
% for ii=1:Nedge
%     if cell_lr(ii,1) == ctgt || cell_lr(ii,2) == ctgt
%         % edge(ii,:)
%         % cell_lr(ii,:)
%         patch('vertices',vtx,'faces',edge(ii,:),'EdgeAlpha',1.0,'Marker','.','Edgecolor','g');
%     end
% end 


% vtxtest = load('io/vtxtest');
% % plot(vtxtest(:,1),vtxtest(:,2),'r.')
% patch('vertices',vtxtest,'faces',connectivity,'EdgeAlpha',0.5,'Marker','.','EdgeColor',[1 0.1 0.1],'MarkerEdgeColor','r');


% vtxtest = load('io/vtxtest1');
% patch('vertices',vtxtest,'faces',connectivity,'EdgeAlpha',0.5,'Marker','.','EdgeColor',[1 0 0]);

% growdir = load('io/growdir1');
% quiver(vtxtest(:,1),vtxtest(:,2),growdir(:,1),growdir(:,2),0.05,'b','maxheadsize',0.01)
% 
% 
% vtxtest = load('io/vtxtest2');
% patch('vertices',vtxtest,'faces',connectivity,'EdgeAlpha',0.5,'Marker','.','EdgeColor',[1 0 0]);
% 
% vtxtest = load('io/vtxtest3');
% patch('vertices',vtxtest,'faces',connectivity,'EdgeAlpha',0.5,'Marker','.','EdgeColor',[1 0 0]);
% 
% vtxtest = load('io/vtxtest4');
% patch('vertices',vtxtest,'faces',connectivity,'EdgeAlpha',0.5,'Marker','.','EdgeColor',[1 0 0]);
% 
% vtxtest = load('io/vtxtest6');
% patch('vertices',vtxtest,'faces',connectivity,'EdgeAlpha',0.5,'Marker','.','EdgeColor',[1 0 0]);
% 
% vtxtest = load('io/vtxtest8');
% patch('vertices',vtxtest,'faces',connectivity,'EdgeAlpha',0.5,'Marker','.','EdgeColor',[1 0 0]);


% 
% vtxtest = load('io/vtxtest');
% patch('vertices',vtxtest,'faces',connectivity,'EdgeAlpha',0.5,'Marker','.','EdgeColor',[0 0 0],'MarkerEdgeColor','r');



% vtxtest = load('io/vtxtest30');
% patch('vertices',vtxtest,'faces',connectivity,'EdgeAlpha',0.5,'Marker','.','EdgeColor',[0 1 0]);
% 
% vtxtest = load('io/vtxtest31');
% patch('vertices',vtxtest,'faces',connectivity,'EdgeAlpha',0.5,'Marker','.','EdgeColor',[0 1 1]);
% 


% vtxtest1 = load('io/vtx_level');
% patch('vertices',vtxtest1,'faces',connectivity,'EdgeAlpha',0.5,'Marker','.','EdgeColor',[1 0 0]);



% vtxtest1 = load('io/vtxtestL1');
% patch('vertices',vtxtest1,'faces',connectivity,'EdgeAlpha',0.5,'Marker','.','EdgeColor',[1 0 0]);

% ndir = load('io/ndir');
% quiver(vtxtest1(:,1),vtxtest1(:,2),ndir(:,1),ndir(:,2),0.01,'b','maxheadsize',0.01)


% vtxtest2 = load('io/vtxtestL2');
% patch('vertices',vtxtest2,'faces',connectivity,'EdgeAlpha',0.5,'Marker','.','EdgeColor',[1 0 0]);
% 
% vtxtest3 = load('io/vtxtestL3');
% patch('vertices',vtxtest3,'faces',connectivity,'EdgeAlpha',0.5,'Marker','.','EdgeColor',[1 0 0]);
% 
% vtxtest4 = load('io/vtxtestL4');
% patch('vertices',vtxtest4,'faces',connectivity,'EdgeAlpha',0.5,'Marker','.','EdgeColor',[1 0 0]);



% for vtgt = 80:110
%     plot([vertices(vtgt,1) vtxtest1(vtgt,1) vtxtest2(vtgt,1) vtxtest3(vtgt,1) vtxtest4(vtgt,1)],[vertices(vtgt,2) vtxtest1(vtgt,2) vtxtest2(vtgt,2) vtxtest3(vtgt,2) vtxtest4(vtgt,2)],'r');
% end

% vtxfull = load('io/vtxfull');
% plot(vtxfull(:,1),vtxfull(:,2),'r*');

% vtxsurf = load('io/vtxsurf');
% plot(vtxsurf(:,1),vtxsurf(:,2),'bo');

% plot(vtxfull(2484,1),vtxfull(2484,2),'b*');
% plot(vtxfull(2612,1),vtxfull(2612,2),'g*');


% vtxtest = load('io/vtxtest');
% plot(vtxtest(:,1),vtxtest(:,2),'r.');
% patch('vertices',vtxtest,'faces',connectivity,'EdgeAlpha',0.5,'Marker','.','EdgeColor',[1 0 0]);
% 
% ndir = load('io/ndir');
% quiver(vtxtest(:,1),vtxtest(:,2),ndir(:,1),ndir(:,2),0.05,'b','maxheadsize',0.01)

% evenbvtx = load('io/evenbvtx');
% plot(evenbvtx(:,1),evenbvtx(:,2),'g*');

% vtgt = 97;
% plot(vertices(vtgt,1),vertices(vtgt,2),'g*');
% vtgt = 193;
% plot(vertices(vtgt,1),vertices(vtgt,2),'r*');
% 
% vtgt = 192;
% plot(vertices(vtgt,1),vertices(vtgt,2),'r*');


%Format
axis equal
axis tight
xlabel('x')
ylabel('y')
hold off

% axis([-0.15    1.15   -0.5    0.5]);

% axis([-0.1699    1.1425   -0.6533    0.6591]);

% axis([0.9823    1.0247    0.0183    0.0608]);

% axis([0.9785    1.0298    0.0332    0.0846]);


% axis([-0.0842    0.0486   -0.0688    0.0640]);

% axis([0.9507    1.0971   -0.0331    0.1136]);

% axis([0.9660    1.0660   -0.0111    0.0891]);

% axis([0.1665    0.2349    0.1439    0.2123]);

% axis([-0.0543    0.0789   -0.0619    0.0714]);

% axis([0.9883    1.0234   -0.0196    0.0155]);
% axis([1.0089    1.0109   -0.0011    0.0009]);

% axis([0.9480    1.0689   -0.0577    0.0635]);
% axis([0.9959    1.0139   -0.0088    0.0092]);


% axis([1.0000    1.0018   -0.0009    0.0009])

% axis([0.9980    1.0043   -0.0031    0.0032]);


% axis([0.9592    1.0344   -0.0394    0.0359]);

% axis([0.8671    1.1031   -0.1265    0.1096]);


% axis([-0.0625    0.3516   -0.1613    0.1573]);







%% Write geometry as loop 

% %Read geometry file 
% [~,~,vertices,connectivity] = import_cell_mesh2d_surface('io/cell_mesh2d_surface.dat');
% Nfcs = length(connectivity(:,1));
% Nvtx = length(vertices(:,1));
% 
% %Find base vertex
% [mval,vbase] = min(vertices(:,1));
% 
% %Build v2f 
% V2F = zeros(Nvtx,2,'int32');
% for ff=1:Nfcs
%     V2F(connectivity(ff,1),2) = ff;
%     V2F(connectivity(ff,2),1) = ff;
% end
% 
% %Accumulate vertices 
% vins = 1;
% vertices_list = zeros(Nvtx,2);
% vertices_list(1,:) = vertices(vbase,:);
% for ii=1:Nvtx-1
%     fnext = V2F(vbase,1);
%     if connectivity(fnext,1) == vbase
%         vnext = connectivity(fnext,2);
%     else
%         vnext = connectivity(fnext,1);
%     end
%     vins = vins + 1;
%     vertices_list(vins,:) = vertices(vnext,:);
%     vbase = vnext;
% end
% 
% %Flip 
% vertices_list(:,1) = 1 - vertices_list(:,1);
% 
% % %Check
% % cla reset
% % hold on
% % plot(vertices_list(:,1),vertices_list(:,2),'r')
% % vtgt = 10;
% % plot(vertices_list(vtgt,1),vertices_list(vtgt,2),'b*')
% % hold off
% 
% %Export
% fid = fopen('io\surface.dat','w+');
% fprintf(fid,'%d \n',Nvtx);
% for ii=1:Nvtx
%     fprintf(fid,'%f %f \n',vertices_list(ii,:));
% end
% fclose(fid);
