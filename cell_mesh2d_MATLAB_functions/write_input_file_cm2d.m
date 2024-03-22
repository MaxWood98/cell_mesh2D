%Function to write cell_mesh2d input file
%Max Wood - mw16116@bristol.ac.uk
%Univeristy of Bristol - Department of Aerospace Engineering

%Version 6.4
%Updated 22-03-2024

%Function -----------------------------------------------------------------
function [] = write_input_file_cm2d(cm2dopt)
fid = fopen(cm2dopt.options_filepath,'w+');
    fprintf(fid,'%s \n','#cell_mesh2d options file (version 0.9.7)');
    fprintf(fid,'%s \n',' ');
    
    fprintf(fid,'%s \n','#=== General Options =================================');
    fprintf(fid,'%s \n','#Toggle console display (yes | no)');
    fprintf(fid,'%s %s \n','condisp =',cm2dopt.condisp);
    fprintf(fid,'%s \n',' ');

    fprintf(fid,'%s \n','#=== Mesh Format Options =============================');
    fprintf(fid,'%s \n','#Type of mesh (cutcell | su2_cutcell | su2_dual | minD_O)');
    fprintf(fid,'%s %s \n','meshtype =',cm2dopt.meshtype);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Mesh inside or outside of geometry');
    fprintf(fid,'%s %s \n','meshinout =',cm2dopt.meshinout);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Surface normal direction switch in to / out of the mesh domain');
    fprintf(fid,'%s %s \n','surfnormdir =',cm2dopt.surface_dir);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Boundary normal direction switch in to / out of the mesh domain');
    fprintf(fid,'%s %s \n','bndrynormdir =',cm2dopt.boundary_dir);
    fprintf(fid,'%s \n',' '); 

    fprintf(fid,'%s \n','#=== Quadtree Options ================================');
    fprintf(fid,'%s \n','#Maximum refinement level');
    fprintf(fid,'%s %d \n','nqtrefine =',cm2dopt.nrefine);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Maximum additional refinement levels in high curvature regions');
    fprintf(fid,'%s %d \n','nboostqtrefine =',cm2dopt.nrefineB);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Maximum number of cells');
    fprintf(fid,'%s %d \n','ncellmax =',cm2dopt.ncell_max);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Refinement adjacency flooding iterations at the first refinement');
    fprintf(fid,'%s %d \n','nadjfloodi =',cm2dopt.nrflood_i);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Refinement adjacency flooding iterations at the final refinement');
    fprintf(fid,'%s %d \n','nadjfloodf =',cm2dopt.nrflood_f);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Refinement adjacency flooding iterations on boosted refinement');
    fprintf(fid,'%s %d \n','nadjfloodb =',cm2dopt.nrflood_b);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Far field distance from object centre');
    fprintf(fid,'%s %f \n','farfielddist =',cm2dopt.fbound);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Object/mesh centre offset (x / y)');
    fprintf(fid,'%s %f \n','offsett_x =',cm2dopt.coffset(1));
    fprintf(fid,'%s %f \n','offsett_y =',cm2dopt.coffset(2));
    fprintf(fid,'%s \n',' '); 

    fprintf(fid,'%s \n','#=== Custom Domain Bound Specification  ==============');
    fprintf(fid,'%s \n','#Set custom mesh domain bounds (yes | no)');
    fprintf(fid,'%s %s \n','forcebounds =',cm2dopt.set_mbounds);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Mesh domain bounds (xmin xmax ymin ymax)');
    fprintf(fid,'%s %E \n','bound_xmin =',cm2dopt.xmin);
    fprintf(fid,'%s %E \n','bound_xmax =',cm2dopt.xmax);
    fprintf(fid,'%s %E \n','bound_ymin =',cm2dopt.ymin);
    fprintf(fid,'%s %E \n','bound_ymax =',cm2dopt.ymax);
    fprintf(fid,'%s \n',' ');

    fprintf(fid,'%s \n','#=== Mesh Cleaning Options ===========================');
    fprintf(fid,'%s \n','#Minimum edge length as a fraction of an undeformed cell edge length at each refienemnt level');
    fprintf(fid,'%s %E \n','eminlength =',cm2dopt.eminlen);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Minimum surface edge length as a fraction of an undeformed cell edge length at each refienemnt level');
    fprintf(fid,'%s %E \n','srfeminlength =',cm2dopt.srfeminlen);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Volume fraction of an undeformed cell at each refinement level below which a cell is classed as a sliver cell');
    fprintf(fid,'%s %s \n','cminvol =',cm2dopt.cminvol);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Allow cleaning of the input surface geometry ');
    fprintf(fid,'%s %s \n','srfinclean =',cm2dopt.srfinclean);
    fprintf(fid,'%s \n',' '); 

    fprintf(fid,'%s \n','#=== Mesh Geometry Intersection Options ==============');
    fprintf(fid,'%s \n','#Maximum number of mesh-geometry intersections on each mesh edge');
    fprintf(fid,'%s %d \n','enintmax =',cm2dopt.enintmax);
    fprintf(fid,'%s \n',' ');     
    fprintf(fid,'%s \n','#Edge-geometry intersection search zone padding as a fraction of the edge length');
    fprintf(fid,'%s %E \n','eintpad =',cm2dopt.eintpad);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Intersection co-incidence tollerance as fraction of item length');
    fprintf(fid,'%s %E \n','intcointol =',cm2dopt.int_coin_tol);
    fprintf(fid,'%s \n',' ');
        
    fprintf(fid,'%s \n','#=== Surface Format Options ==========================');
    fprintf(fid,'%s \n','#Geometry surface type (simplified | exact)');
    fprintf(fid,'%s %s \n','surftype =',cm2dopt.surftype);
    fprintf(fid,'%s \n',' ');  
    fprintf(fid,'%s \n','#Surface curvature multiplier');
    fprintf(fid,'%s %f \n','scurvmult =',cm2dopt.surfRcurvM);
    fprintf(fid,'%s \n',' ');  
    fprintf(fid,'%s \n','#Number of vertices used to estimate local surface curvaturer');
    fprintf(fid,'%s %d \n','scurvnpnt =',cm2dopt.surfRcurvNpts);
    fprintf(fid,'%s \n',' ');  

    fprintf(fid,'%s \n','#=== Inflation Layer Options =========================');
    fprintf(fid,'%s \n','#Toggle construction of an inflation layer (yes | no)');
    fprintf(fid,'%s %s \n','build_inflayer =',cm2dopt.build_inflayer);
    fprintf(fid,'%s \n',' ');  
    fprintf(fid,'%s \n','#Inflation layer approximate height');
    fprintf(fid,'%s %f \n','inflayer_height =',cm2dopt.inflayer_height);
    fprintf(fid,'%s \n',' ');  
    fprintf(fid,'%s \n','#Number of layers within the inflation layer (determine automatically if zero)');
    fprintf(fid,'%s %d \n','inflayer_nlayer =',cm2dopt.inflayer_nlayer);
    fprintf(fid,'%s \n',' ');    
    fprintf(fid,'%s \n','#Inflation layer differencing weight (0->1)');
    fprintf(fid,'%s %f \n','inflayer_wd =',cm2dopt.inflayer_wd);
    fprintf(fid,'%s \n',' ');  
    fprintf(fid,'%s \n','#Inflation layer convex differencing penalty (0->1)');
    fprintf(fid,'%s %f \n','inflayer_cvxdp =',cm2dopt.inflayer_cvxdp);
    fprintf(fid,'%s \n',' ');  
    fprintf(fid,'%s \n','#Inflation layer evening weight (0->1)');
    fprintf(fid,'%s %f \n','inflayer_we =',cm2dopt.inflayer_we);
    fprintf(fid,'%s \n',' ');  
    fprintf(fid,'%s \n','#Inflation layer evening bias base value');
    fprintf(fid,'%s %f \n','inflayer_ebcbase =',cm2dopt.inflayer_ebcbase);
    fprintf(fid,'%s \n',' ');  
    fprintf(fid,'%s \n','#Inflation layer evening number of flood iterations');
    fprintf(fid,'%s %d \n','inflayer_enflood =',cm2dopt.inflayer_enflood);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Inflation layer evening number of sub-iterations');
    fprintf(fid,'%s %d \n','inflayer_ensubiter =',cm2dopt.inflayer_ensubiter);
    fprintf(fid,'%s \n',' ');  
    fprintf(fid,'%s \n','#Inflation layer evening length ratio bound (>=1)');
    fprintf(fid,'%s %f \n','inflayer_lreb =',cm2dopt.inflayer_lreb);
    fprintf(fid,'%s \n',' ');  

    fprintf(fid,'%s \n','#=== Mesh Smoothing Options ==========================');
    fprintf(fid,'%s \n','#Near surface smoothing type (none | laplacian)');
    fprintf(fid,'%s %s \n','smoothingtype =',cm2dopt.nsstype);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Vertex selection flooding iterations from surfaces');
    fprintf(fid,'%s %d \n','smoothingnflood = ',cm2dopt.nsvtxflood);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Smoothing iterations');
    fprintf(fid,'%s %d \n','smoothingniter =',cm2dopt.nlpsmooth);
    fprintf(fid,'%s \n',' '); 
    
    fprintf(fid,'%s \n','#=== ADtree Options ==================================');
    fprintf(fid,'%s \n','#Padding size of adtree search bounding boxes as a multiple of cell edge length');
    fprintf(fid,'%s %f \n','adtpadding =',cm2dopt.adtree_spad);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#AD tree maximum depth as the number of dimension cycles (tree is 4d)');
    fprintf(fid,'%s %d \n','adtndimcycle =',cm2dopt.adtree_maxd);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#AD tree minimum divisible node size');
    fprintf(fid,'%s %d \n','adtmindivnsize =',cm2dopt.adtree_mindivnsize);
    fprintf(fid,'%s \n',' '); 
        
    fprintf(fid,'%s \n','#=== Gradient Projection Options =====================');
    fprintf(fid,'%s \n','#Construct and export volume to surface gradient interpolation (yes | no)');
    fprintf(fid,'%s %s \n','glinkconstructexp =',cm2dopt.glink_con);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Gradient linking type (rbf or int)');
    fprintf(fid,'%s %s \n','glinktype =',cm2dopt.glinktype);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Number of nearest neighbours to use for RBF volume to surface gradient interpolation');
    fprintf(fid,'%s %d \n','glinkrbfnpnt =',cm2dopt.glink_nnn);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Number of vertices each side used to smooth the gradient at each surface vertex');
    fprintf(fid,'%s %d \n','glinknpntsmooth =',cm2dopt.glink_nsmooth);
    fprintf(fid,'%s \n',' ');

    fprintf(fid,'%s \n','#=== Radial Basis Function Options ===================');
    fprintf(fid,'%s \n','#RBF interpolation support radius as a multiple of the maximum distance between points');
    fprintf(fid,'%s %f \n','rbfrsup =',cm2dopt.RBF_rsup);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#RBF interpolation relaxation distance as fraction of the maximum distance between points');
    fprintf(fid,'%s %f \n','rbfrrelax =',cm2dopt.RBF_relaxD);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#RBF interpolation relaxation parameter');
    fprintf(fid,'%s %f \n','rbfrelaxp =',cm2dopt.RBF_relaxP);
    fprintf(fid,'%s \n',' ');

    fprintf(fid,'%s \n','#=== Boundary Condition Options ======================');
    fprintf(fid,'%s \n','#Set custom boundary conditions in specifed regions (yes | no)');
    fprintf(fid,'%s %s \n','setcustombcs =',cm2dopt.set_custom_bc);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Remove any region of the mesh connected to a far field boundary condition (yes | no)');
    fprintf(fid,'%s %s \n','remffczones =',cm2dopt.rem_ffzones);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Remove any region of the mesh connected to a non custom set boundary condition (yes | no)');
    fprintf(fid,'%s %s \n','remncczones =',cm2dopt.rem_nczones);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Remove any isolated region of the mesh connected only to a wall boundary condition (yes | no)');
    fprintf(fid,'%s %s \n','remwallonlyzones =',cm2dopt.rem_iszones);
    fprintf(fid,'%s \n',' '); 
fclose(fid);
end