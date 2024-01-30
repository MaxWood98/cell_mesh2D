%Function to write cell_mesh2d input file (v2)
%Max Wood - mw16116@bristol.ac.uk
%Univeristy of Bristol - Department of Aerospace Engineering

%Version 3.0
%Updated 07-11-2023

%Function -----------------------------------------------------------------
function [] = write_input_file_cm2d(cm2dop)
fid = fopen('io\cell_mesh2d_options.dat','w+');
    fprintf(fid,'%s \n','#cell_mesh2d V2 options file (version 0.7.0)');
    fprintf(fid,'%s \n',' ');
    
    fprintf(fid,'%s \n','#=== General Options =================================');
    fprintf(fid,'%s \n','#Toggle console display (1 = yes || 0 = no)');
    fprintf(fid,'%d \n',cm2dop.condisp);
    fprintf(fid,'%s \n',' ');

    fprintf(fid,'%s \n','#=== Mesh Format Options =============================');
    fprintf(fid,'%s \n','#Type of mesh (0 = cutcell | 1 = minD block mesh)');
    fprintf(fid,'%d \n',cm2dop.meshtype);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Mesh inside or outside of geometry');
    fprintf(fid,'%s \n',cm2dop.meshinout);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Surface normal direction switch in to / out of the mesh domain');
    fprintf(fid,'%s \n',cm2dop.surface_dir);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Boundary normal direction switch in to / out of the mesh domain');
    fprintf(fid,'%s \n',cm2dop.boundary_dir);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Mesh output format (cutcell / su2_cutcell / su2_dual)');
    fprintf(fid,'%s \n',cm2dop.meshfrmat);
    fprintf(fid,'%s \n',' '); 

    fprintf(fid,'%s \n','#=== Quadtree Options ================================');
    fprintf(fid,'%s \n','#Maximum refinement level');
    fprintf(fid,'%d \n',cm2dop.nrefine);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Maximum additional refinement levels in high curvature regions');
    fprintf(fid,'%d \n',cm2dop.nrefineB);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Maximum number of cells');
    fprintf(fid,'%d \n',cm2dop.ncell_max);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Refinement adjacency flooding iterations at the first refinement');
    fprintf(fid,'%d \n',cm2dop.nrflood_i);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Refinement adjacency flooding iterations at the final refinement');
    fprintf(fid,'%d \n',cm2dop.nrflood_f);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Refinement adjacency flooding iterations on boosted refinement');
    fprintf(fid,'%d \n',cm2dop.nrflood_b);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Far field distance from object centre');
    fprintf(fid,'%f \n',cm2dop.fbound);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Object/mesh centre offset (x / y)');
    fprintf(fid,'%f \n',cm2dop.coffset(1));
    fprintf(fid,'%f \n',cm2dop.coffset(2));
    fprintf(fid,'%s \n',' '); 

    fprintf(fid,'%s \n','#=== Custom Domain Bound Specification  ==============');
    fprintf(fid,'%s \n','#Enable custom mesh domain bounds (1 = yes | 0 = no)');
    fprintf(fid,'%d \n',cm2dop.set_mbounds);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Mesh domain bounds (xmin xmax ymin ymax)');
    fprintf(fid,'%E %E %E %E \n',cm2dop.xmin,cm2dop.xmax,cm2dop.ymin,cm2dop.ymax);
    fprintf(fid,'%s \n',' ');

    fprintf(fid,'%s \n','#=== Mesh Cleaning Options ===========================');
    fprintf(fid,'%s \n','#Minimum edge length as a fraction of an undeformed cell edge length at each refienemnt level');
    fprintf(fid,'%E \n',cm2dop.eminlen);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Volume fraction of an undeformed cell at each refinement level below which a cell is classed as a sliver cell');
    fprintf(fid,'%E \n',cm2dop.cminvol);
    fprintf(fid,'%s \n',' '); 

    fprintf(fid,'%s \n','#=== Mesh Geometry Intersection Options ==============');
    fprintf(fid,'%s \n','#Maximum number of mesh-geometry intersections on each mesh edge');
    fprintf(fid,'%d \n',cm2dop.enintmax);
    fprintf(fid,'%s \n',' ');     
    fprintf(fid,'%s \n','#Edge-geometry intersection search zone padding as a fraction of the edge length');
    fprintf(fid,'%E \n',cm2dop.eintpad);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Intersection co-incidence tollerance as fraction of item length');
    fprintf(fid,'%E \n',cm2dop.int_coin_tol);
    fprintf(fid,'%s \n',' ');
        
    fprintf(fid,'%s \n','#=== Surface Format Options ==========================');
    fprintf(fid,'%s \n','#Geometry surface type (0 = simplified | 1 = exact)');
    fprintf(fid,'%d \n',cm2dop.surftype);
    fprintf(fid,'%s \n',' ');  
    fprintf(fid,'%s \n','#Surface curvature multiplier');
    fprintf(fid,'%f \n',cm2dop.surfRcurvM);
    fprintf(fid,'%s \n',' ');  
    fprintf(fid,'%s \n','#Number of vertices used to estimate local surface curvaturer');
    fprintf(fid,'%d \n',cm2dop.surfRcurvNpts);
    fprintf(fid,'%s \n',' ');  

    fprintf(fid,'%s \n','#=== Mesh Smoothing Options ==========================');
    fprintf(fid,'%s \n','#Near surface smoothing type (0 = none | 1 = Laplacian)');
    fprintf(fid,'%d \n',cm2dop.nsstype);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Vertex selection flooding iterations from surfaces');
    fprintf(fid,'%d \n',cm2dop.nsvtxflood);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Smoothing iterations');
    fprintf(fid,'%d \n',cm2dop.nlpsmooth);
    fprintf(fid,'%s \n',' '); 
    
    fprintf(fid,'%s \n','#=== ADtree Options ==================================');
    fprintf(fid,'%s \n','#Padding size of adtree search bounding boxes as a multiple of cell edge length');
    fprintf(fid,'%f \n',cm2dop.adtree_spad);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#AD tree maximum depth as the number of dimension cycles (tree is 4d)');
    fprintf(fid,'%d \n',cm2dop.adtree_maxd);
    fprintf(fid,'%s \n',' '); 
    
    fprintf(fid,'%s \n','#=== Gradient Interpolation Options ==================');
    fprintf(fid,'%s \n','#Construct volume to surface gradient interpolation (1 = yes | 0 = no)');
    fprintf(fid,'%d \n',cm2dop.glink_con);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Gradient linking type (rbf or int)');
    fprintf(fid,'%s \n',cm2dop.glinktype);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Number of nearest neighbours to use for volume to surface gradient interpolation');
    fprintf(fid,'%d \n',cm2dop.glink_nnn);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Number of vertices each side used to smooth the gradient at each surface vertex');
    fprintf(fid,'%d \n',cm2dop.glink_nsmooth);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#RBF interpolation smoothing relaxation parameter');
    fprintf(fid,'%f \n',cm2dop.glink_RBF_relax);
    fprintf(fid,'%s \n',' ');

    fprintf(fid,'%s \n','#=== Boundary Condition Options ======================');
    fprintf(fid,'%s \n','#Set custom boundary conditions in specifed regions (1 = yes | 0 = no)');
    fprintf(fid,'%d \n',cm2dop.set_custom_bc);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Remove any region of the mesh connected to a far field boundary condition');
    fprintf(fid,'%d \n',cm2dop.rem_ffzones);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Remove any region of the mesh connected to a non custom set boundary condition');
    fprintf(fid,'%d \n',cm2dop.rem_nczones);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Remove any isolated region of the mesh connected only to a wall boundary condition');
    fprintf(fid,'%d \n',cm2dop.rem_iszones);
    fprintf(fid,'%s \n',' '); 
fclose(fid);
end