!cell_mesh2d gradient coupling module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 0.8
!Updated 24-07-2023

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
integer(in) :: ee,ii,jj,vv,nn,kk,ff,pp,aa
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

            !Cases
            if (vsins .LT. cm2dopt%glink_nnn) then !If too few point have been found 

                !Reset arrays 
                R(:,:) = 0.0d0 
                Ri(:,:) = 0.0d0
                surf2volRBF(:) = 0.0d0 

                !Build local dependance matrix and set support radius 
                do ii=1,vsins
                    v1 = vsurf_gnn(ii)
                    do jj=1,vsins
                        v2 = vsurf_gnn(jj)
                        R(ii,jj) = norm2(volume_mesh%vertices(v2,:) - volume_mesh%vertices(v1,:)) 
                    end do 
                end do 
                Rs = 50.0d0*maxval(R(1:vsins,1:vsins))
                do ii=1,vsins
                    do jj=1,vsins
                        R(ii,jj) = wendlandc2(R(ii,jj),Rs)
                    end do 
                end do 

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
                    do ii=1,cm2dopt%glink_nnn
                        v1 = vsurf_gnn(ii)
                        do jj=1,cm2dopt%glink_nnn
                            v2 = vsurf_gnn(jj)
                            R(ii,jj) = norm2(volume_mesh%vertices(v2,:) - volume_mesh%vertices(v1,:)) 
                        end do 
                    end do 
                    Rs = 50.0d0*maxval(R)
                    do ii=1,cm2dopt%glink_nnn
                        do jj=1,cm2dopt%glink_nnn
                            R(ii,jj) = wendlandc2(R(ii,jj),Rs)
                        end do 
                    end do 

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


end module cellmesh2d_gradient_coupling_mod