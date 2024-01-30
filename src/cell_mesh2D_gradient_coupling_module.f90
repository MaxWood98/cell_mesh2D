!cell_mesh2d gradient coupling module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 0.9
!Updated 07-11-2023

!Module
module cellmesh2d_gradient_coupling_mod
use  cellmesh2d_surface_mod
use cellmesh2d_geometry_mod
contains 

!Gradient coupling matrix construction subroutine ===========================
subroutine construct_surfvol_grad_coupling(volume_mesh,surface_mesh,Ndim,node_minDIVsize,global_target_pad,cm2dopt)
implicit none 

!Variables - Import
integer(in) :: Ndim,node_minDIVsize
real(dp) :: global_target_pad
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh
type(cm2d_options) :: cm2dopt

!Variables - Local 
integer(in) :: ee,ii,vv,nn,kk,ff,pp,aa
integer(in) :: Nesurf,nselected,eminD,v1,v2,nvsurf_sel,nvsurf_selN,nvbase,etgt,vtgt,nadd,vsins,vselect,vb,va,Npts
integer(in) :: vsurf_select(volume_mesh%nvtx),vsurf_tag(volume_mesh%nvtx),vsurf_base(volume_mesh%nvtx)
integer(in) :: vsurf_dtagged(volume_mesh%nvtx),vsurf_gnn(cm2dopt%glink_nnn),smoothinterp(2*cm2dopt%glink_nsmooth+1)
integer(in) :: edge_surf(volume_mesh%nedge,2),v2e_vsurf(volume_mesh%nvtx,2)
integer(in), dimension(:), allocatable :: node_select
real(dp) :: lpad,dist_c,dist_n,vs_maxd,Rs,zxmin,zxmax,zymin,zymax,zzmin,zzmax,s_base
real(dp) :: vl1(2),vl2(2),vid(3),vsurf_dist(volume_mesh%nvtx),smooth_s(2*cm2dopt%glink_nsmooth+1)
real(dp) :: surf2volRBF(cm2dopt%glink_nnn)
real(dp) :: R(cm2dopt%glink_nnn,cm2dopt%glink_nnn),Ri(cm2dopt%glink_nnn,cm2dopt%glink_nnn)
real(dp), dimension(:,:), allocatable :: tvtx
type(tree_data) :: sv_adtree

!Allocate surface interpolation structure 
allocate(volume_mesh%vsinterp(surface_mesh%nvtx))

!Extract surface mesh within volume mesh 
Nesurf = 0 
edge_surf(:,:) = 0
do ee=1,volume_mesh%nedge
    if ((volume_mesh%edge(ee,3) == -1) .OR. (volume_mesh%edge(ee,4) == -1)) then 
        Nesurf = Nesurf + 1
        edge_surf(Nesurf,:) = volume_mesh%edge(ee,1:2)
    end if 
end do 

!Build v2e for this surface sub-mesh
v2e_vsurf(:,:) = 0 
do ee=1,Nesurf
    v1 = edge_surf(ee,1)
    v2 = edge_surf(ee,2)
    v2e_vsurf(v1,2) = ee 
    v2e_vsurf(v2,1) = ee 
end do 

!Build adtree on this volume surface mesh 
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '    {constructing AD-tree on surface within volume mesh}'
end if
allocate(tvtx(Nesurf,4))
do ii=1,Nesurf
    tvtx(ii,1) = minval(volume_mesh%vertices(edge_surf(ii,:),1)) !xmin
    tvtx(ii,2) = minval(volume_mesh%vertices(edge_surf(ii,:),2)) !ymin
    tvtx(ii,3) = maxval(volume_mesh%vertices(edge_surf(ii,:),1)) !xmax
    tvtx(ii,4) = maxval(volume_mesh%vertices(edge_surf(ii,:),2)) !ymax
end do
call build_ADtree(sv_adtree,ndim,cm2dopt%ADTmax_depth,node_minDIVsize,tvtx,global_target_pad,cm2dopt%dispt)
nselected = 0
allocate(node_select(sv_adtree%nnode))
node_select(:) = 0 

!Build local interpolation for each surface point
zxmin = 0.0d0 
zxmax = 0.0d0 
zymin = 0.0d0 
zymax = 0.0d0 
zzmin = 0.0d0 
zzmax = 0.0d0
R(:,:) = 0.0d0
Ri(:,:) = 0.0d0
vsurf_tag(:) = 0 
vsurf_base(:) = 0 
vsurf_select(:) = 0 
vsurf_dtagged(:) = 0
vsurf_dist(:) = 0.0d0  
do vv=1,surface_mesh%nvtx
    if (surface_mesh%vtx_active(vv) == 1) then 

        !Set search length as 4x the maximum cell size on this vertex 
        lpad = 4.0d0*surface_mesh%vtx_rsearch(vv)

        !Find the closest volume mesh surface edge to the current surface mesh vertex
        zxmin = surface_mesh%vertices(vv,1) - lpad !tgt bounding box -> xmin
        zxmax = surface_mesh%vertices(vv,1) + lpad !tgt bounding box -> xmax
        zymin = surface_mesh%vertices(vv,2) - lpad !tgt bounding box -> ymin
        zymax = surface_mesh%vertices(vv,2) + lpad !tgt bounding box -> ymax
        call search_ADtree(nselected,node_select,sv_adtree,zxmin,zxmax,zymin,zymax,zzmin,zzmax)
        eminD = 0 
        dist_c = 2.0d0*cm2dopt%far_field_bound
        do nn=1,nselected
            do kk=1,sv_adtree%tree(node_select(nn))%nentry
        
                !Verticies of the ends of the surface segment in the volume mesh
                vl1(:) = volume_mesh%vertices(edge_surf(sv_adtree%tree(node_select(nn))%entry(kk),1),:)
                vl2(:) = volume_mesh%vertices(edge_surf(sv_adtree%tree(node_select(nn))%entry(kk),2),:)

                !Test distance from surface point to this segment 
                vid = min_dist_point_to_edge(vl1,vl2,surface_mesh%vertices(vv,:))
                dist_n = vid(3)
                if (dist_n .LE. dist_c) then 
                    eminD = sv_adtree%tree(node_select(nn))%entry(kk)
                    dist_c = dist_n
                end if 
            end do 
        end do 
        !print *, vv,eminD,dist_c

        !If non zero target edge (else if zero then this edge is outside the mesh domain so can be ignored)
        if (eminD .NE. 0) then 

            !Find glink_nnn nearest neighbours superset
            nvbase = 2
            nvsurf_sel = 2 
            vsurf_base(1) = edge_surf(eminD,1)
            vsurf_base(2) = edge_surf(eminD,2)
            vsurf_select(1) = edge_surf(eminD,1)
            vsurf_select(2) = edge_surf(eminD,2)
            vsurf_dist(1) = norm2(volume_mesh%vertices(vsurf_select(1),:) - surface_mesh%vertices(vv,:))
            vsurf_dist(2) = norm2(volume_mesh%vertices(vsurf_select(2),:) - surface_mesh%vertices(vv,:))
            vs_maxd = maxval(vsurf_dist(1:2))
            vsurf_tag(vsurf_select(1:2)) = 1
            do ff=1,volume_mesh%nvtx !Flood to add valid adjacent new points 
                nadd = 0 
                nvsurf_selN = nvsurf_sel
                do pp=1,nvsurf_sel

                    !Edge 1
                    etgt = v2e_vsurf(vsurf_select(pp),1)
                    if (etgt .GT. 0) then 
                        if (edge_surf(etgt,1) == vsurf_select(pp)) then 
                            vtgt = edge_surf(etgt,2)
                        else
                            vtgt = edge_surf(etgt,1)
                        end if 
                        if (vsurf_tag(vtgt) == 0) then !New vertex so add if valid
                            if (nvsurf_sel .LE. cm2dopt%glink_nnn) then !If not enough points yet found 
                                nadd = nadd + 1
                                nvbase = nvbase + 1
                                vsurf_base(nvbase) = vtgt
                                vsurf_tag(vtgt) = 1
                                nvsurf_selN = nvsurf_selN + 1
                                vsurf_select(nvsurf_selN) = vtgt 
                                vsurf_dist(nvsurf_selN) = norm2(volume_mesh%vertices(vtgt,:) - surface_mesh%vertices(vv,:))
                                if (vsurf_dist(nvsurf_selN) .GT. vs_maxd) then !Update maximum distance if required
                                    vs_maxd = vsurf_dist(nvsurf_selN)
                                end if
                            elseif (norm2(volume_mesh%vertices(vtgt,:) - surface_mesh%vertices(vv,:)) .LE. vs_maxd) then !If enough points have been found but this point is closer than the current maximum distance point
                                nadd = nadd + 1
                                nvbase = nvbase + 1
                                vsurf_base(nvbase) = vtgt
                                vsurf_tag(vtgt) = 1
                                nvsurf_selN = nvsurf_selN + 1
                                vsurf_select(nvsurf_selN) = vtgt 
                                vsurf_dist(nvsurf_selN) = norm2(volume_mesh%vertices(vtgt,:) - surface_mesh%vertices(vv,:))
                                if (vsurf_dist(nvsurf_selN) .GT. vs_maxd) then !Update maximum distance if required
                                    vs_maxd = vsurf_dist(nvsurf_selN)
                                end if
                            end if 
                        end if 
                    end if 

                    !Edge 2
                    etgt = v2e_vsurf(vsurf_select(pp),2)
                    if (etgt .GT. 0) then 
                        if (edge_surf(etgt,1) == vsurf_select(pp)) then 
                            vtgt = edge_surf(etgt,2)
                        else
                            vtgt = edge_surf(etgt,1)
                        end if 
                        if (vsurf_tag(vtgt) == 0) then !New vertex so add if valid
                            if (nvsurf_sel .LE. cm2dopt%glink_nnn) then !If not enough points yet found 
                                nadd = nadd + 1
                                nvbase = nvbase + 1
                                vsurf_base(nvbase) = vtgt
                                vsurf_tag(vtgt) = 1
                                nvsurf_selN = nvsurf_selN + 1
                                vsurf_select(nvsurf_selN) = vtgt 
                                vsurf_dist(nvsurf_selN) = norm2(volume_mesh%vertices(vtgt,:) - surface_mesh%vertices(vv,:))
                                if (vsurf_dist(nvsurf_selN) .GT. vs_maxd) then !Update maximum distance if required
                                    vs_maxd = vsurf_dist(nvsurf_selN)
                                end if
                            elseif (norm2(volume_mesh%vertices(vtgt,:) - surface_mesh%vertices(vv,:)) .LE. vs_maxd) then !If enough points have been found but this point is closer than the current maximum distance point
                                nadd = nadd + 1
                                nvbase = nvbase + 1
                                vsurf_base(nvbase) = vtgt
                                vsurf_tag(vtgt) = 1
                                nvsurf_selN = nvsurf_selN + 1
                                vsurf_select(nvsurf_selN) = vtgt 
                                vsurf_dist(nvsurf_selN) = norm2(volume_mesh%vertices(vtgt,:) - surface_mesh%vertices(vv,:))
                                if (vsurf_dist(nvsurf_selN) .GT. vs_maxd) then !Update maximum distance if required
                                    vs_maxd = vsurf_dist(nvsurf_selN)
                                end if
                            end if 
                        end if 
                    end if 
                end do 

                !Set new point count 
                nvsurf_sel = nvsurf_selN

                !Exit if no new points added 
                if (nadd == 0) then 
                    exit 
                end if 
            end do 
            vsurf_tag(vsurf_select(1:nvsurf_sel)) = 0 !Reset vertex tags 

            !Find actual glink_nnn nearest neighbours
            vsins = 0 
            vsurf_gnn(:) = 0 
            vsurf_dtagged(1:nvsurf_sel) = 0 
            do pp=1,nvsurf_sel
                vselect = minloc(vsurf_dist(1:nvsurf_sel),1,mask = vsurf_dtagged(1:nvsurf_sel) == 0)
                vsurf_dtagged(vselect) = 1
                vsurf_dist(vselect) = 2.0d0*cm2dopt%far_field_bound
                vtgt = vsurf_select(vselect)
                vsins = vsins + 1
                vsurf_gnn(vsins) = vtgt 
                if (vsins == cm2dopt%glink_nnn) then !Exit when all points are found 
                    exit 
                end if 
            end do 

            !Reset arrays 
            Ri(:,:) = 0.0d0
            surf2volRBF(:) = 0.0d0 

            !Cases
            if (vsins .LT. cm2dopt%glink_nnn) then !If too few point have been found 
                
                !Build local dependance matrix and set support radius 
                call build_RBF_influence(R,Rs,vsins,vsurf_gnn,volume_mesh%vertices,cm2dopt)

                !Invert dependance matrix 
                call matinv(R(1:vsins,1:vsins),Ri(1:vsins,1:vsins),vsins)

                !Build RBF dependance of surface point on all selected volume points 
                do ii=1,vsins
                    dist_c = norm2(volume_mesh%vertices(vsurf_gnn(ii),:) - surface_mesh%vertices(vv,:))
                    surf2volRBF(ii) = wendlandc2(dist_c,Rs)
                end do 

                !Set number of points 
                Npts = vsins
            else !If full number of points found
                if (cm2dopt%glink_nnn .GT. 1) then !More than one surface point requested

                    !Build local dependance matrix and set support radius 
                    call build_RBF_influence(R,Rs,cm2dopt%glink_nnn,vsurf_gnn,volume_mesh%vertices,cm2dopt)

                    !Invert dependance matrix 
                    call matinv(R,Ri,cm2dopt%glink_nnn)

                    !Build RBF dependance of surface point on all selected volume points 
                    do ii=1,cm2dopt%glink_nnn
                        dist_c = norm2(volume_mesh%vertices(vsurf_gnn(ii),:) - surface_mesh%vertices(vv,:))
                        surf2volRBF(ii) = wendlandc2(dist_c,Rs)
                    end do 

                    !Set number of points 
                    Npts = cm2dopt%glink_nnn 
                else !One surface point requested
                    surf2volRBF(1) = 1.0d0 
                    Ri(1,1) = 1.0d0 
                    Npts = 1
                end if 
            end if 

            !Store interpolation structure 
            volume_mesh%vsinterp(vv)%npnts_vi = Npts
            allocate(volume_mesh%vsinterp(vv)%vol_pnts(cm2dopt%glink_nnn))
            allocate(volume_mesh%vsinterp(vv)%surf2volRBF(cm2dopt%glink_nnn))
            allocate(volume_mesh%vsinterp(vv)%Ri(cm2dopt%glink_nnn,cm2dopt%glink_nnn))
            volume_mesh%vsinterp(vv)%vol_pnts(:) = 0.0d0 
            volume_mesh%vsinterp(vv)%surf2volRBF(:) = 0.0d0 
            volume_mesh%vsinterp(vv)%Ri(:,:) = 0.0d0 
            volume_mesh%vsinterp(vv)%vol_pnts(1:Npts) = vsurf_gnn(1:Npts)
            volume_mesh%vsinterp(vv)%surf2volRBF(1:Npts) = surf2volRBF(1:Npts)
            volume_mesh%vsinterp(vv)%Ri(1:Npts,1:Npts) = Ri(1:Npts,1:Npts)
        else
            volume_mesh%vsinterp(vv)%npnts_vi = 0 
            allocate(volume_mesh%vsinterp(vv)%vol_pnts(cm2dopt%glink_nnn))
            allocate(volume_mesh%vsinterp(vv)%surf2volRBF(cm2dopt%glink_nnn))
            allocate(volume_mesh%vsinterp(vv)%Ri(cm2dopt%glink_nnn,cm2dopt%glink_nnn))
            volume_mesh%vsinterp(vv)%vol_pnts(:) = 0
            volume_mesh%vsinterp(vv)%surf2volRBF(:) = 0
            volume_mesh%vsinterp(vv)%Ri(:,:) = 0.0d0 
        end if 
    else
        volume_mesh%vsinterp(vv)%npnts_vi = 0 
        allocate(volume_mesh%vsinterp(vv)%vol_pnts(cm2dopt%glink_nnn))
        allocate(volume_mesh%vsinterp(vv)%surf2volRBF(cm2dopt%glink_nnn))
        allocate(volume_mesh%vsinterp(vv)%Ri(cm2dopt%glink_nnn,cm2dopt%glink_nnn))
        volume_mesh%vsinterp(vv)%vol_pnts(:) = 0
        volume_mesh%vsinterp(vv)%surf2volRBF(:) = 0
        volume_mesh%vsinterp(vv)%Ri(:,:) = 0.0d0 
    end if 
end do 

!Add linked vertices to smooth each surface vertex 
do vv=1,surface_mesh%nvtx
    if (surface_mesh%vtx_active(vv) == 1) then 

        !Build linked interpolation set 
        smoothinterp(:) = 0 
        vb = vv 
        do aa=1,cm2dopt%glink_nsmooth !previous vertices
            va = get_previous_vertex(surface_mesh,vb)
            smoothinterp(cm2dopt%glink_nsmooth+1-aa) = va
            vb = va 
        end do  
        smoothinterp(cm2dopt%glink_nsmooth+1) = vv 
        vb = vv 
        do aa=1,cm2dopt%glink_nsmooth !next vertices
            va = get_next_vertex(surface_mesh,vb)
            smoothinterp(cm2dopt%glink_nsmooth+1+aa) = va
            vb = va 
        end do 

        !Find surface coordinate of each interpolation point and the base point 
        smooth_s(:) = 0.0d0 
        do aa=2,2*cm2dopt%glink_nsmooth+1
            smooth_s(aa) = smooth_s(aa-1) + norm2(surface_mesh%vertices(smoothinterp(aa),:) - &
            surface_mesh%vertices(smoothinterp(aa-1),:))
        end do 
        s_base = smooth_s(cm2dopt%glink_nsmooth+1)

        !Allocate and store smoothing points
        allocate(volume_mesh%vsinterp(vv)%surf_smooth_pnts(2*cm2dopt%glink_nsmooth+1))
        allocate(volume_mesh%vsinterp(vv)%surf_smoothRBF(2*cm2dopt%glink_nsmooth+1))
        volume_mesh%vsinterp(vv)%surf_smooth_pnts(:) = smoothinterp(:)

        !Evaluate and store RBF weighting for each point along the surface
        if (cm2dopt%glink_nsmooth .GE. 1) then 
            Rs = 1.1d0*smooth_s(2*cm2dopt%glink_nsmooth+1)
            do aa=1,2*cm2dopt%glink_nsmooth+1
                volume_mesh%vsinterp(vv)%surf_smoothRBF(aa) = wendlandc2(abs(smooth_s(aa)-s_base),Rs)
            end do 
        else
            volume_mesh%vsinterp(vv)%surf_smoothRBF(1) = 1.0d0 
        end if 
    else    

        !Allocate and set to zero smoothing points
        allocate(volume_mesh%vsinterp(vv)%surf_smooth_pnts(2*cm2dopt%glink_nsmooth+1))
        allocate(volume_mesh%vsinterp(vv)%surf_smoothRBF(2*cm2dopt%glink_nsmooth+1))
        volume_mesh%vsinterp(vv)%surf_smooth_pnts(:) = 0
        volume_mesh%vsinterp(vv)%surf_smoothRBF(:) = 0 
    end if 
end do 
return 
end subroutine construct_surfvol_grad_coupling




!Gradient projection subroutine (RBF interpolation) ===========================
subroutine project_gradients_RBF(gradient_surf,gradient_vol,volume_mesh,surface_mesh,Ndim,node_minDIVsize,global_target_pad,cm2dopt)
implicit none 

!Variables - Import
integer(in) :: Ndim,node_minDIVsize
real(dp) :: global_target_pad
real(dp), dimension(:,:) :: gradient_vol
real(dp), dimension(:,:), allocatable :: gradient_surf
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh
type(cm2d_options) :: cm2dopt

!Variables - Local 
integer(in) :: ee,ii,vv,nn,kk,ff,pp,aa
integer(in) :: Nesurf,nselected,eminD,v1,v2,nvsurf_sel,nvsurf_selN,nvbase,ftgt,etgt,vtgt,nadd,vsins,vselect,vb,va,Npts
integer(in) :: vsurf_select(volume_mesh%nvtx),vsurf_tag(volume_mesh%nvtx),vsurf_base(volume_mesh%nvtx)
integer(in) :: vsurf_dtagged(volume_mesh%nvtx),vsurf_gnn(cm2dopt%glink_nnn),smoothinterp(2*cm2dopt%glink_nsmooth+1)
integer(in) :: edge_surf(volume_mesh%nedge,2),v2e_vsurf(volume_mesh%nvtx,2)
integer(in), dimension(:), allocatable :: node_select
real(dp) :: lpad,dist_c,dist_n,vs_maxd,Rs,zxmin,zxmax,zymin,zymax,zzmin,zzmax,s_base
real(dp) :: vl1(2),vl2(2),vid(3),vsurf_dist(volume_mesh%nvtx),smooth_s(2*cm2dopt%glink_nsmooth+1)
real(dp) :: surf2volRBF(cm2dopt%glink_nnn)
real(dp) :: R(cm2dopt%glink_nnn,cm2dopt%glink_nnn),Ri(cm2dopt%glink_nnn,cm2dopt%glink_nnn),surf_smoothRBF(2*cm2dopt%glink_nsmooth+1)
real(dp), dimension(:,:), allocatable :: tvtx,gradient_surf_temp
type(tree_data) :: sv_adtree

!Allocate and initialise surface gradient structures 
allocate(gradient_surf(surface_mesh%nvtx,2))
allocate(gradient_surf_temp(surface_mesh%nvtx,2))
gradient_surf(:,:) = 0.0d0 
gradient_surf_temp(:,:) = 0.0d0 

!Allocate surface interpolation structure 
allocate(volume_mesh%vsinterp(surface_mesh%nvtx))

!Extract surface mesh within volume mesh 
Nesurf = 0 
edge_surf(:,:) = 0
do ee=1,volume_mesh%nedge
    if ((volume_mesh%edge(ee,3) == -1) .OR. (volume_mesh%edge(ee,4) == -1)) then 
        Nesurf = Nesurf + 1
        edge_surf(Nesurf,:) = volume_mesh%edge(ee,1:2)
    end if 
end do 

!Build v2e for this surface sub-mesh
v2e_vsurf(:,:) = 0 
do ee=1,Nesurf
    v1 = edge_surf(ee,1)
    v2 = edge_surf(ee,2)
    v2e_vsurf(v1,2) = ee 
    v2e_vsurf(v2,1) = ee 
end do 

!Build adtree on this volume surface mesh 
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '    {constructing AD-tree on surface within volume mesh}'
end if
allocate(tvtx(Nesurf,4))
do ii=1,Nesurf
    tvtx(ii,1) = minval(volume_mesh%vertices(edge_surf(ii,:),1)) !xmin
    tvtx(ii,2) = minval(volume_mesh%vertices(edge_surf(ii,:),2)) !ymin
    tvtx(ii,3) = maxval(volume_mesh%vertices(edge_surf(ii,:),1)) !xmax
    tvtx(ii,4) = maxval(volume_mesh%vertices(edge_surf(ii,:),2)) !ymax
end do
call build_ADtree(sv_adtree,ndim,cm2dopt%ADTmax_depth,node_minDIVsize,tvtx,global_target_pad,cm2dopt%dispt)
nselected = 0
allocate(node_select(sv_adtree%nnode))
node_select(:) = 0 

!Build local interpolation for each surface point
zxmin = 0.0d0 
zxmax = 0.0d0 
zymin = 0.0d0 
zymax = 0.0d0 
zzmin = 0.0d0 
zzmax = 0.0d0
R(:,:) = 0.0d0
Ri(:,:) = 0.0d0
vsurf_tag(:) = 0 
vsurf_base(:) = 0 
vsurf_select(:) = 0 
vsurf_dtagged(:) = 0
vsurf_dist(:) = 0.0d0  
do vv=1,surface_mesh%nvtx

    !Find the distance to the farthest connected surface vertex from this vertex
    dist_c = 0.0d0 
    do aa=1,2
        ftgt =  surface_mesh%v2f(vv,aa)
        dist_n = norm2(surface_mesh%vertices(surface_mesh%faces(ftgt,2),:) - &
        surface_mesh%vertices(surface_mesh%faces(ftgt,1),:))
        if (dist_n .GT. dist_c) then 
            dist_c = dist_n
        end if 
    end do 

    !Set search length as 4x the distance to the farthest other surface vertex
    lpad = 4.0d0*dist_c

    !Find the closest volume mesh surface edge to the current surface mesh vertex
    zxmin = surface_mesh%vertices(vv,1) - lpad !tgt bounding box -> xmin
    zxmax = surface_mesh%vertices(vv,1) + lpad !tgt bounding box -> xmax
    zymin = surface_mesh%vertices(vv,2) - lpad !tgt bounding box -> ymin
    zymax = surface_mesh%vertices(vv,2) + lpad !tgt bounding box -> ymax
    call search_ADtree(nselected,node_select,sv_adtree,zxmin,zxmax,zymin,zymax,zzmin,zzmax)
    eminD = 0 
    dist_c = 2.0d0*cm2dopt%far_field_bound
    do nn=1,nselected
        do kk=1,sv_adtree%tree(node_select(nn))%nentry
    
            !Verticies of the ends of the surface segment in the volume mesh
            vl1(:) = volume_mesh%vertices(edge_surf(sv_adtree%tree(node_select(nn))%entry(kk),1),:)
            vl2(:) = volume_mesh%vertices(edge_surf(sv_adtree%tree(node_select(nn))%entry(kk),2),:)

            !Test distance from surface point to this segment 
            vid = min_dist_point_to_edge(vl1,vl2,surface_mesh%vertices(vv,:))
            dist_n = vid(3)
            if (dist_n .LE. dist_c) then 
                eminD = sv_adtree%tree(node_select(nn))%entry(kk)
                dist_c = dist_n
            end if 
        end do 
    end do 
    ! print *, vv,eminD,dist_c,lpad

    !If non zero target edge (else if zero then this edge is outside the mesh domain so can be ignored)
    if (eminD .NE. 0) then 

        !Find glink_nnn nearest neighbours superset
        nvbase = 2
        nvsurf_sel = 2 
        vsurf_base(1) = edge_surf(eminD,1)
        vsurf_base(2) = edge_surf(eminD,2)
        vsurf_select(1) = edge_surf(eminD,1)
        vsurf_select(2) = edge_surf(eminD,2)
        vsurf_dist(1) = norm2(volume_mesh%vertices(vsurf_select(1),:) - surface_mesh%vertices(vv,:))
        vsurf_dist(2) = norm2(volume_mesh%vertices(vsurf_select(2),:) - surface_mesh%vertices(vv,:))
        vs_maxd = maxval(vsurf_dist(1:2))
        vsurf_tag(vsurf_select(1:2)) = 1
        do ff=1,volume_mesh%nvtx !Flood to add valid adjacent new points 
            nadd = 0 
            nvsurf_selN = nvsurf_sel
            do pp=1,nvsurf_sel

                !Edge 1
                etgt = v2e_vsurf(vsurf_select(pp),1)
                if (etgt .GT. 0) then 
                    if (edge_surf(etgt,1) == vsurf_select(pp)) then 
                        vtgt = edge_surf(etgt,2)
                    else
                        vtgt = edge_surf(etgt,1)
                    end if 
                    if (vsurf_tag(vtgt) == 0) then !New vertex so add if valid
                        if (nvsurf_sel .LE. cm2dopt%glink_nnn) then !If not enough points yet found 
                            nadd = nadd + 1
                            nvbase = nvbase + 1
                            vsurf_base(nvbase) = vtgt
                            vsurf_tag(vtgt) = 1
                            nvsurf_selN = nvsurf_selN + 1
                            vsurf_select(nvsurf_selN) = vtgt 
                            vsurf_dist(nvsurf_selN) = norm2(volume_mesh%vertices(vtgt,:) - surface_mesh%vertices(vv,:))
                            if (vsurf_dist(nvsurf_selN) .GT. vs_maxd) then !Update maximum distance if required
                                vs_maxd = vsurf_dist(nvsurf_selN)
                            end if
                        elseif (norm2(volume_mesh%vertices(vtgt,:) - surface_mesh%vertices(vv,:)) .LE. vs_maxd) then !If enough points have been found but this point is closer than the current maximum distance point
                            nadd = nadd + 1
                            nvbase = nvbase + 1
                            vsurf_base(nvbase) = vtgt
                            vsurf_tag(vtgt) = 1
                            nvsurf_selN = nvsurf_selN + 1
                            vsurf_select(nvsurf_selN) = vtgt 
                            vsurf_dist(nvsurf_selN) = norm2(volume_mesh%vertices(vtgt,:) - surface_mesh%vertices(vv,:))
                            if (vsurf_dist(nvsurf_selN) .GT. vs_maxd) then !Update maximum distance if required
                                vs_maxd = vsurf_dist(nvsurf_selN)
                            end if
                        end if 
                    end if 
                end if 

                !Edge 2
                etgt = v2e_vsurf(vsurf_select(pp),2)
                if (etgt .GT. 0) then 
                    if (edge_surf(etgt,1) == vsurf_select(pp)) then 
                        vtgt = edge_surf(etgt,2)
                    else
                        vtgt = edge_surf(etgt,1)
                    end if 
                    if (vsurf_tag(vtgt) == 0) then !New vertex so add if valid
                        if (nvsurf_sel .LE. cm2dopt%glink_nnn) then !If not enough points yet found 
                            nadd = nadd + 1
                            nvbase = nvbase + 1
                            vsurf_base(nvbase) = vtgt
                            vsurf_tag(vtgt) = 1
                            nvsurf_selN = nvsurf_selN + 1
                            vsurf_select(nvsurf_selN) = vtgt 
                            vsurf_dist(nvsurf_selN) = norm2(volume_mesh%vertices(vtgt,:) - surface_mesh%vertices(vv,:))
                            if (vsurf_dist(nvsurf_selN) .GT. vs_maxd) then !Update maximum distance if required
                                vs_maxd = vsurf_dist(nvsurf_selN)
                            end if
                        elseif (norm2(volume_mesh%vertices(vtgt,:) - surface_mesh%vertices(vv,:)) .LE. vs_maxd) then !If enough points have been found but this point is closer than the current maximum distance point
                            nadd = nadd + 1
                            nvbase = nvbase + 1
                            vsurf_base(nvbase) = vtgt
                            vsurf_tag(vtgt) = 1
                            nvsurf_selN = nvsurf_selN + 1
                            vsurf_select(nvsurf_selN) = vtgt 
                            vsurf_dist(nvsurf_selN) = norm2(volume_mesh%vertices(vtgt,:) - surface_mesh%vertices(vv,:))
                            if (vsurf_dist(nvsurf_selN) .GT. vs_maxd) then !Update maximum distance if required
                                vs_maxd = vsurf_dist(nvsurf_selN)
                            end if
                        end if 
                    end if 
                end if 
            end do 

            !Set new point count 
            nvsurf_sel = nvsurf_selN

            !Exit if no new points added 
            if (nadd == 0) then 
                exit 
            end if 
        end do 
        vsurf_tag(vsurf_select(1:nvsurf_sel)) = 0 !Reset vertex tags 

        !Find actual glink_nnn nearest neighbours
        vsins = 0 
        vsurf_gnn(:) = 0 
        vsurf_dtagged(1:nvsurf_sel) = 0 
        do pp=1,nvsurf_sel
            vselect = minloc(vsurf_dist(1:nvsurf_sel),1,mask = vsurf_dtagged(1:nvsurf_sel) == 0)
            vsurf_dtagged(vselect) = 1
            vsurf_dist(vselect) = 2.0d0*cm2dopt%far_field_bound
            vtgt = vsurf_select(vselect)
            vsins = vsins + 1
            vsurf_gnn(vsins) = vtgt 
            if (vsins == cm2dopt%glink_nnn) then !Exit when all points are found 
                exit 
            end if 
        end do 

        !Reset arrays 
        Ri(:,:) = 0.0d0
        surf2volRBF(:) = 0.0d0 

        !Cases
        if (vsins .LT. cm2dopt%glink_nnn) then !If too few point have been found 

            !Build local dependance matrix and set support radius 
            call build_RBF_influence(R,Rs,vsins,vsurf_gnn,volume_mesh%vertices,cm2dopt)

            !Invert dependance matrix 
            call matinv(R(1:vsins,1:vsins),Ri(1:vsins,1:vsins),vsins)

            !Build RBF dependance of surface point on all selected volume points 
            do ii=1,vsins
                dist_c = norm2(volume_mesh%vertices(vsurf_gnn(ii),:) - surface_mesh%vertices(vv,:))
                surf2volRBF(ii) = wendlandc2(dist_c,Rs)
            end do 

            !Set number of points 
            Npts = vsins
        else !If full number of points found
            if (cm2dopt%glink_nnn .GT. 1) then !More than one surface point requested

                !Build local dependance matrix and set support radius 
                call build_RBF_influence(R,Rs,cm2dopt%glink_nnn,vsurf_gnn,volume_mesh%vertices,cm2dopt)

                !Invert dependance matrix 
                call matinv(R,Ri,cm2dopt%glink_nnn)

                !Build RBF dependance of surface point on all selected volume points 
                do ii=1,cm2dopt%glink_nnn
                    dist_c = norm2(volume_mesh%vertices(vsurf_gnn(ii),:) - surface_mesh%vertices(vv,:))
                    surf2volRBF(ii) = wendlandc2(dist_c,Rs)
                end do 

                !Set number of points 
                Npts = cm2dopt%glink_nnn 
            else !One surface point requested
                surf2volRBF(1) = 1.0d0 
                Ri(1,1) = 1.0d0 
                Npts = 1
            end if 
        end if 

        !Project to this surface vertex 
        gradient_surf(vv,:) = project_vtx_vgrad_2_surfgrad(gradient_vol,surf2volRBF(1:Npts),Ri,vsurf_gnn(1:Npts),Npts)
        volume_mesh%vsinterp(vv)%npnts_vi = Npts
    else !no face found so ignore projection
        !do nothing 
    end if 
end do 

!Return if no smoothing 
if (cm2dopt%glink_nsmooth == 0) then 
    return 
end if

!Store surface gradients in the temporary structure
gradient_surf_temp(:,:) = gradient_surf(:,:)

!Display
if (cm2dopt%dispt == 1) then
    write(*,'(A)') '    {smoothing surface gradients}'
end if

!Add linked vertices to smooth each surface vertex 
do vv=1,surface_mesh%nvtx
    if (volume_mesh%vsinterp(vv)%npnts_vi .GT. 0) then !smooth this vertex

        !Build linked interpolation set 
        smoothinterp(:) = 0 
        vb = vv 
        do aa=1,cm2dopt%glink_nsmooth !previous vertices
            va = get_previous_vertex(surface_mesh,vb)
            smoothinterp(cm2dopt%glink_nsmooth+1-aa) = va
            vb = va 
        end do  
        smoothinterp(cm2dopt%glink_nsmooth+1) = vv 
        vb = vv 
        do aa=1,cm2dopt%glink_nsmooth !next vertices
            va = get_next_vertex(surface_mesh,vb)
            smoothinterp(cm2dopt%glink_nsmooth+1+aa) = va
            vb = va 
        end do 

        !Find surface coordinate of each interpolation point and the base point 
        smooth_s(:) = 0.0d0 
        do aa=2,2*cm2dopt%glink_nsmooth+1
            smooth_s(aa) = smooth_s(aa-1) + norm2(surface_mesh%vertices(smoothinterp(aa),:) - &
            surface_mesh%vertices(smoothinterp(aa-1),:))
        end do 
        s_base = smooth_s(cm2dopt%glink_nsmooth+1)

        !Build RBF values at smoothing points 
        surf_smoothRBF(:) = 0.0d0 
        if (cm2dopt%glink_nsmooth+1 .GE. 1) then 
            Rs = 1.1d0*maxval(smooth_s(:))
            do aa=1,2*cm2dopt%glink_nsmooth+1
                surf_smoothRBF(aa) = wendlandc2(abs(smooth_s(aa)-s_base),Rs)
            end do 
        else
            surf_smoothRBF(1) = 1.0d0 
        end if 

        !Apply smoothing to this vertex 
        gradient_surf(vv,:) = smooth_surface_gradient(gradient_surf_temp,surf_smoothRBF,smoothinterp,2*cm2dopt%glink_nsmooth+1)
    else !no interpolation here so ignore smoothing 
        !do nothing 
    end if 
end do 
return 
end subroutine project_gradients_RBF




!Gradient projection subroutine (direct edge interpolation) ===========================
subroutine project_gradients_INT(gradient_surf,gradient_vol,volume_mesh,surface_mesh,cm2dopt)
implicit none 

!Variables - Import
real(dp), dimension(:,:) :: gradient_vol
real(dp), dimension(:,:), allocatable :: gradient_surf
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh
type(cm2d_options) :: cm2dopt

!Variables - Local
integer(in) :: ii
integer(in) :: etgt,vtgt,ev1,ev2
real(dp) :: sedge_frac

!Allocate and initialise surface gradient structure
allocate(gradient_surf(surface_mesh%nvtx,2))
gradient_surf(:,:) = 0.0d0 

!Accumulate projected gradients to each surface vertex
do ii=1,volume_mesh%nvtx_surf 

    !Target volume mesh vertex, surface edge and surface edge fraction
    vtgt = volume_mesh%surf_vtx(ii)
    etgt = volume_mesh%surf_vtx_seg(ii)
    sedge_frac = volume_mesh%surf_vtx_segfrac(ii)

    !Surface vertices on this edge
    ev1 = surface_mesh%faces(etgt,1)
    ev2 = surface_mesh%faces(etgt,2)

    !Accumulate gradient at rate 1-sedge_frac to ev1 and sedge_frac to ev2
    gradient_surf(ev1,1) = gradient_surf(ev1,1) + (1.0d0 - sedge_frac)*gradient_vol(vtgt,1)
    gradient_surf(ev1,2) = gradient_surf(ev1,2) + (1.0d0 - sedge_frac)*gradient_vol(vtgt,2)
    gradient_surf(ev2,1) = gradient_surf(ev2,1) + sedge_frac*gradient_vol(vtgt,1)
    gradient_surf(ev2,2) = gradient_surf(ev2,2) + sedge_frac*gradient_vol(vtgt,2)
end do 
return 
end subroutine project_gradients_INT




!Subroutine to build RBF influence ===========================
subroutine build_RBF_influence(R,R_sup,Npoint,point_list,vertices,cm2dopt)
implicit none 

!Variables - Import
integer(in) :: Npoint 
integer(in), dimension(:) :: point_list
real(dp) :: R_sup
real(dp), dimension(:,:) :: R,vertices 
type(cm2d_options) :: cm2dopt

!Variables - Local
integer(in) :: ii,jj 
integer(in) :: v1,v2
real(dp) :: mindist(Npoint)

!Populate distance matrix
R(:,:) = 0.0d0 
mindist(:) = ieee_value(1.0d0,IEEE_POSITIVE_INF)
do ii=1,Npoint
    v1 = point_list(ii)
    do jj=1,Npoint
        v2 = point_list(jj)
        R(ii,jj) = norm2(vertices(v2,:) - vertices(v1,:)) 
        if ((R(ii,jj) .GT. 0.0d0) .AND. (R(ii,jj) .LT. mindist(ii))) then 
            mindist(ii) = R(ii,jj)
        end if 
    end do 
end do 

!Set support radius 
R_sup = 50.0d0*maxval(R(1:Npoint,1:Npoint))

!Evaluate RBF values for each distance entry
do ii=1,Npoint
    do jj=1,Npoint
        R(ii,jj) = wendlandc2(R(ii,jj),R_sup)
    end do
end do 

!Relax interpolation at interior points 
if ((Npoint .GT. 2) .AND. (cm2dopt%RBF_relax .NE. 0.0d0)) then 
    do ii=2,Npoint-1
        R(ii,ii) = R(ii,ii) + cm2dopt%RBF_relax/(mindist(ii)*R_sup)
    end do 
end if 
return 
end subroutine build_RBF_influence




!Smooth gradient at vertex function ===========================
function smooth_surface_gradient(gradient_surf_unsmooth,smoothRBF,surfpoints,Nsurfpoint) result(sgrad_smooth)
implicit none 

!Result
real(dp) :: sgrad_smooth(2)

!Variables - Import
integer(in) :: Nsurfpoint
integer(in), dimension(:) :: surfpoints
real(dp), dimension(:) :: smoothRBF
real(dp), dimension(:,:) :: gradient_surf_unsmooth

!Variables - Local
integer(in) :: ii 
real(dp) :: weight_sum

!Find smoothed gradient at this vertex
weight_sum = 0.0d0 
sgrad_smooth(:) = 0.0d0 
do ii=1,Nsurfpoint
    weight_sum = weight_sum + smoothRBF(ii)
    sgrad_smooth(:) = sgrad_smooth(:) + gradient_surf_unsmooth(surfpoints(ii),:)*smoothRBF(ii)
end do 
sgrad_smooth(:) = sgrad_smooth(:)/weight_sum
return 
end function smooth_surface_gradient




!Project gradient to vertex function ===========================
function project_vtx_vgrad_2_surfgrad(gradient_vol,surf2volRBF,Ri,volpoints,Nvolp) result(vsgrad)
implicit none 

!Result
real(dp) :: vsgrad(2)

!Variables - Import
integer(in) :: Nvolp
integer(in), dimension(:) :: volpoints
real(dp), dimension(:) :: surf2volRBF
real(dp), dimension(:,:) :: Ri,gradient_vol

!Variables - Local
integer(in) :: ii 
real(dp) :: gamma_xy(Nvolp,2),gradient_vol_loc(Nvolp,2)

!Populate local gradient vector
do ii=1,Nvolp
    gradient_vol_loc(ii,:) = gradient_vol(volpoints(ii),:)
end do 

!Solve gamma coefficients 
gamma_xy(:,1) = matmul(Ri(1:Nvolp,1:Nvolp),gradient_vol_loc(:,1))
gamma_xy(:,2) = matmul(Ri(1:Nvolp,1:Nvolp),gradient_vol_loc(:,2))

!Construct interpolated gradient 
vsgrad(1) = dot_product(gamma_xy(:,1),surf2volRBF(:))
vsgrad(2) = dot_product(gamma_xy(:,2),surf2volRBF(:))
return 
end function project_vtx_vgrad_2_surfgrad


end module cellmesh2d_gradient_coupling_mod