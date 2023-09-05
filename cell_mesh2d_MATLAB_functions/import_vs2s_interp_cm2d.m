%Function to read cell_mesh2d volume to surface interpolation file 
%Max Wood - mw16116@bristol.ac.uk
%Univeristy of Bristol - Department of Aerospace Engineering

%Version 2.1
%Updated 24-07-2023

%Function -----------------------------------------------------------------
function [vs2s_interp] = import_vs2s_interp_cm2d(filename)
    fid = fopen(filename,'r');
    idata = textscan(fid,'%d %d %d',1);
    Nvsurf = idata{1,1};
    Npinterp = idata{1,2};
    Npsmooth = idata{1,3};
    vs2s_interp.Nvsurf = Nvsurf;
    vs2s_interp.Npinterp = Npinterp;
    vs2s_interp.Npsmooth = Npsmooth;
    vs2s_interp.idata = cell(Nvsurf,1);
    vpfrmat = '%d';
    addvpfrmt = ' %d';
    for ii=1:Npinterp-1
        vpfrmat = [vpfrmat addvpfrmt];
    end
    spfrmat = '%d';
    addspfrmt = ' %d';
    for ii=1:Npsmooth-1
        spfrmat = [spfrmat addspfrmt];
    end
    
    smatfrmat = '%f';
    addmatfrmt = ' %f';
    for ii=1:Npsmooth-1
        smatfrmat = [smatfrmat addmatfrmt];
    end

    matfrmat = '%f';
    addmatfrmt = ' %f';
    for ii=1:Npinterp-1
        matfrmat = [matfrmat addmatfrmt];
    end
    for vv=1:Nvsurf
        %textscan(fid,'%s',1) %skip line
        spidx = textscan(fid,'%d',1); %vertex index
        Npnts = textscan(fid,'%d',1); %number of interpolation points
        vs2s_interp.idata{vv}.Npnts = Npnts{1,1};
        smoothpoints = textscan(fid,spfrmat,1);
        vs2s_interp.idata{vv}.smoothpoints = zeros(Npsmooth,1);
        for ii=1:Npsmooth
            vs2s_interp.idata{vv}.smoothpoints(ii) = smoothpoints{ii};
        end
        smoothpointsRBF = textscan(fid,smatfrmat,1);
        vs2s_interp.idata{vv}.smoothpointsRBF = zeros(Npsmooth,1);
        for ii=1:Npsmooth
            vs2s_interp.idata{vv}.smoothpointsRBF(ii) = smoothpointsRBF{ii};
        end
        volpoints = textscan(fid,vpfrmat,1);
        vs2s_interp.idata{vv}.volpoints = zeros(Npinterp,1);
        for ii=1:Npinterp
            vs2s_interp.idata{vv}.volpoints(ii) = volpoints{ii};
        end
        vs2s_interp.idata{vv}.surf2volRBF = zeros(Npinterp,1);
        surf2volRBF = textscan(fid,matfrmat,1);
        for ii=1:Npinterp
            vs2s_interp.idata{vv}.surf2volRBF(ii) = surf2volRBF{ii};
        end
        matrix = textscan(fid,matfrmat,Npinterp);
        vs2s_interp.idata{vv}.matrix = zeros(Npinterp,Npinterp);
        for jj=1:Npinterp
            matcol = matrix{1,jj};
            for ii=1:Npinterp
                vs2s_interp.idata{vv}.matrix(ii,jj) = matcol(ii);
            end
        end
    end
    fclose(fid);
end