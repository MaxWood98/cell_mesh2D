!cell_mesh2d distance field module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 0.1
!Updated 07-02-2023

!Module 
module cellmesh2d_distF_mod
use cellmesh2d_geometry_mod
contains 


!Preprocess mesh for distance field construction subroutine ================================================
subroutine preprocess_mesh(meshdat,volume_mesh)
implicit none 

!Variables - Import 
type(meshdata) :: meshdat 
type(vol_mesh_data) :: volume_mesh

!Variables - Local  
integer(in) :: ii,vv,cc,ee,v1,v2,c3,c4,exist,ctgt,e1,e2 
integer(in) :: v2e(volume_mesh%nvtx,4)
real(dp) :: Aedg,Ledge,dx,dy,cosang,twopi
real(dp) :: cell_lwedge(volume_mesh%ncell),vs1(2),vs2(2)

!Define 2pi 
twopi = 8.0d0*atan(1.0d0)

!Allocate
allocate(meshdat%vtx_vlnc(volume_mesh%nvtx))
allocate(meshdat%vtx_v2v(volume_mesh%nvtx,4))
allocate(meshdat%vtx_v2e(volume_mesh%nvtx,4))
allocate(meshdat%edgedx(volume_mesh%nedge))
allocate(meshdat%edgedy(volume_mesh%nedge))
allocate(meshdat%edgedxd(volume_mesh%nedge))
allocate(meshdat%edgedyd(volume_mesh%nedge))
allocate(meshdat%cell_area(volume_mesh%ncell))
allocate(meshdat%cell_midp(volume_mesh%ncell,2))
allocate(meshdat%edge_midp(volume_mesh%nedge,2))
meshdat%vtx_vlnc(:) = 0 
meshdat%vtx_v2v(:,:) = 0 
meshdat%vtx_v2e(:,:) = 0 
meshdat%edgedx(:) = 0.0d0 
meshdat%edgedy(:) = 0.0d0 
meshdat%edgedxd(:) = 0.0d0 
meshdat%edgedyd(:) = 0.0d0 
meshdat%cell_area(:) = 0.0d0 
meshdat%cell_midp(:,:) = 0.0d0 
meshdat%edge_midp(:,:) = 0.0d0 

!Evaluate main parameters
cell_lwedge(:) = 0.0d0 
do ii=1,volume_mesh%nedge
    v1 = volume_mesh%edge(ii,1)
    v2 = volume_mesh%edge(ii,2)
    c3 = volume_mesh%edge(ii,3)
    c4 = volume_mesh%edge(ii,4)
    Aedg = Asegment(volume_mesh%vertices(v1,:),volume_mesh%vertices(v2,:))
    meshdat%edge_midp(ii,1) = 0.5d0*(volume_mesh%vertices(v1,1) + volume_mesh%vertices(v2,1))
    meshdat%edge_midp(ii,2) = 0.5d0*(volume_mesh%vertices(v1,2) + volume_mesh%vertices(v2,2))
    dx = volume_mesh%vertices(v2,1) - volume_mesh%vertices(v1,1)
    dy = volume_mesh%vertices(v2,2) - volume_mesh%vertices(v1,2)
    meshdat%edgedx(ii) = dx 
    meshdat%edgedy(ii) = dy
    Ledge = sqrt(dx**2 + dy**2)
    meshdat%cell_area(c4) = meshdat%cell_area(c4) - Aedg
    meshdat%cell_midp(c4,1) = meshdat%cell_midp(c4,1) + meshdat%edge_midp(ii,1)*Ledge
    meshdat%cell_midp(c4,2) = meshdat%cell_midp(c4,2) + meshdat%edge_midp(ii,2)*Ledge
    cell_lwedge(c4) = cell_lwedge(c4) + Ledge
    if (c3 .GT. 0) then  
        meshdat%cell_area(c3) = meshdat%cell_area(c3) + Aedg
        meshdat%cell_midp(c3,1) = meshdat%cell_midp(c3,1) + meshdat%edge_midp(ii,1)*Ledge
        meshdat%cell_midp(c3,2) = meshdat%cell_midp(c3,2) + meshdat%edge_midp(ii,2)*Ledge
        cell_lwedge(c3) = cell_lwedge(c3) + Ledge
    end if 
end do 
meshdat%cell_midp(:,1) = meshdat%cell_midp(:,1)/cell_lwedge(:)
meshdat%cell_midp(:,2) = meshdat%cell_midp(:,2)/cell_lwedge(:)
do ii=1,volume_mesh%nedge

    !Edge properties
    c3 = volume_mesh%edge(ii,3)
    c4 = volume_mesh%edge(ii,4)

    !Dual edge direction
    if (c3 .GT. 0) then 
        dx = meshdat%cell_midp(c3,1) - meshdat%cell_midp(c4,1)
        dy = meshdat%cell_midp(c3,2) - meshdat%cell_midp(c4,2)
    else
        ! dx = 2.0d0*(meshdat%edge_midp(ii,1) - meshdat%cell_midp(c4,1))
        ! dy = 2.0d0*(meshdat%edge_midp(ii,2) - meshdat%cell_midp(c4,2))
        dx = meshdat%edge_midp(ii,1) - meshdat%cell_midp(c4,1)
        dy = meshdat%edge_midp(ii,2) - meshdat%cell_midp(c4,2)
    end if 
    meshdat%edgedxd(ii) = dx
    meshdat%edgedyd(ii) = dy
end do 

!Find vertex to cell mapping in mesh Ra
allocate(meshdat%vtx_2_cell(volume_mesh%nvtx,4))
meshdat%vtx_2_cell(:,:) = 0 
do ii=1,volume_mesh%nedge

    !Edge adjacent cells 
    c3 = volume_mesh%edge(ii,3)
    c4 = volume_mesh%edge(ii,4)

    !Add cells to vertices on this edge
    do vv=1,2
        v1 = volume_mesh%edge(ii,vv)
        exist = 0 
        do cc=1,4
            if (meshdat%vtx_2_cell(v1,cc) == c4) then 
                exist = 1
            end if
        end do 
        if (exist == 0) then 
            do cc=1,4
                if (meshdat%vtx_2_cell(v1,cc) == 0) then 
                    meshdat%vtx_2_cell(v1,cc) = c4 
                    exit
                end if
            end do 
        end if
        if (c3 .GT. 0) then 

            exist = 0 
            do cc=1,4
                if (meshdat%vtx_2_cell(v1,cc) == c3) then 
                    exist = 1
                end if
            end do 
            if (exist == 0) then 
                do cc=1,4
                    if (meshdat%vtx_2_cell(v1,cc) == 0) then 
                        meshdat%vtx_2_cell(v1,cc) = c3 
                        exit
                    end if
                end do 
            end if
        end if 
    end do 
end do 

!Build v2e for each vertex in Ra
v2e(:,:) = 0 
do ii=1,volume_mesh%nedge

    !Vertices on this edge 
    v1 = volume_mesh%edge(ii,1)
    v2 = volume_mesh%edge(ii,2)

    !Add to v2e
    do vv=1,4
        if (v2e(v1,vv) == 0) then 
            v2e(v1,vv) = ii
            exit
        end if
    end do 
    do vv=1,4
        if (v2e(v2,vv) == 0) then 
            v2e(v2,vv) = ii
            exit
        end if 
    end do 
end do 
meshdat%vtx_v2e(:,:) = v2e(:,:) 

!Build valence and v2v for each vertex 
do ii=1,volume_mesh%nedge

    !Vertices on this edge 
    v1 = volume_mesh%edge(ii,1)
    v2 = volume_mesh%edge(ii,2)

    !Accumulate
    meshdat%vtx_vlnc(v1) = meshdat%vtx_vlnc(v1) + 1
    meshdat%vtx_v2v(v1,meshdat%vtx_vlnc(v1)) = v2 
    meshdat%vtx_vlnc(v2) = meshdat%vtx_vlnc(v2) + 1
    meshdat%vtx_v2v(v2,meshdat%vtx_vlnc(v2)) = v1 
end do

!Assign cell weightings to each vertex
allocate(meshdat%cell_vtx_W(volume_mesh%nvtx,4))
meshdat%cell_vtx_W(:,:) = 0.0d0 
do vv=1,volume_mesh%nvtx

    !Each cell on this vertex
    do cc=1,4
        if (meshdat%vtx_2_cell(vv,cc) .NE. 0) then 

            !Cell 
            ctgt = meshdat%vtx_2_cell(vv,cc)

            !Find edges on this cell 
            e1 = 0
            e2 = 0 
            do ee=1,4
                if (v2e(vv,ee) .NE. 0) then 
                    if ((volume_mesh%edge(v2e(vv,ee),3) == ctgt) .OR.  (volume_mesh%edge(v2e(vv,ee),4) == ctgt)) then 
                        e1 = v2e(vv,ee)
                        exit 
                    end if 
                end if 
            end do 
            do ee=1,4
                if ((v2e(vv,ee) .NE. 0) .AND. (v2e(vv,ee) .NE. e1)) then 
                    if ((volume_mesh%edge(v2e(vv,ee),3) == ctgt) .OR.  (volume_mesh%edge(v2e(vv,ee),4) == ctgt)) then 
                        e2 = v2e(vv,ee)
                        exit 
                    end if 
                end if 
            end do

            !Find span vectors
            if (volume_mesh%edge(e1,1) == vv) then 
                v1 = volume_mesh%edge(e1,2)
            else
                v1 = volume_mesh%edge(e1,1)
            end if
            if (volume_mesh%edge(e2,1) == vv) then 
                v2 = volume_mesh%edge(e2,2)
            else
                v2 = volume_mesh%edge(e2,1)
            end if
            vs1(:) = volume_mesh%vertices(v1,:) - volume_mesh%vertices(vv,:)
            vs2(:) = volume_mesh%vertices(v2,:) - volume_mesh%vertices(vv,:)
            vs1(:) = vs1(:)/sqrt(vs1(1)**2 + vs1(2)**2)
            vs2(:) = vs2(:)/sqrt(vs2(1)**2 + vs2(2)**2)

            !Find weighting 
            cosang = vs1(1)*vs2(1) + vs1(2)*vs2(2)
            meshdat%cell_vtx_W(vv,cc) = abs(acos(cosang)/twopi)
            ! meshdat%cell_vtx_W(vv,cc) = meshdat%cell_area(ctgt)*abs(acos(cosang)/twopi)
        end if 
    end do
end do 
do vv=1,volume_mesh%nvtx
    meshdat%cell_vtx_W(vv,:) = meshdat%cell_vtx_W(vv,:)/sum(meshdat%cell_vtx_W(vv,:))
end do 

!Find cell2edge mapping 
allocate(meshdat%cell_nedge(volume_mesh%ncell))
meshdat%cell_nedge(:) = 0
do ii=1,volume_mesh%nedge 
    c3 = volume_mesh%edge(ii,3)
    c4 = volume_mesh%edge(ii,4)
    if (c3 .GT. 0) then 
        meshdat%cell_nedge(c3) = meshdat%cell_nedge(c3) + 1
    end if
    if (c4 .GT. 0) then 
        meshdat%cell_nedge(c4) = meshdat%cell_nedge(c4) + 1
    end if
end do 
allocate(meshdat%cell2edge(volume_mesh%ncell,maxval(meshdat%cell_nedge(:))))
meshdat%cell_nedge(:) = 0 
meshdat%cell2edge(:,:) = 0 
do ii=1,volume_mesh%nedge 
    c3 = volume_mesh%edge(ii,3)
    c4 = volume_mesh%edge(ii,4)
    if (c3 .GT. 0) then 
        meshdat%cell_nedge(c3) = meshdat%cell_nedge(c3) + 1
        meshdat%cell2edge(c3,meshdat%cell_nedge(c3)) = ii 
    end if
    if (c4 .GT. 0) then 
        meshdat%cell_nedge(c4) = meshdat%cell_nedge(c4) + 1
        meshdat%cell2edge(c4,meshdat%cell_nedge(c4)) = ii 
    end if
end do 

!Identify active boundary conditions 
allocate(meshdat%bc_active(6))
meshdat%bc_active(:) = 0 
do ii=1,volume_mesh%nedge 
    c3 = volume_mesh%edge(ii,3)
    if (c3 .LT. 0) then 
        meshdat%bc_active(abs(c3)) = 1
    end if 
end do 
return 
end subroutine preprocess_mesh




!Build distance field function ================================================
subroutine build_distance_field(cell_d,grad_d,grad2_d,volume_mesh,meshdat,BCactive,cm2dopt)
implicit none 

!System data 
type(meshdata) :: meshdat 
type(cm2d_options) :: cm2dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Import 
integer(in) :: BCactive(volume_mesh%nedge)
real(dp) :: cell_d(volume_mesh%ncell)
real(dp) :: grad_d(volume_mesh%ncell,2)
real(dp) :: grad2_d(volume_mesh%ncell)

!Variables - Local 
integer(in) :: solveiter,ee
integer(in) :: c3,c4,niter
real(dp) :: rc3,rc4,re,edx,edy,gxf,gyf,gr4x,gr4y,gr3x,gr3y,kd,CFL,conresbound
real(dp) :: dfx,dfy,dfm2,dfm,gre_fvx,gre_fvy,flx_x,flx_y,flx_edg,dres
real(dp) :: cell_phi(volume_mesh%ncell)
real(dp) :: phi_grad2(volume_mesh%ncell)
real(dp) :: phi_grad_mag(volume_mesh%ncell)
real(dp) :: phi_grad(volume_mesh%ncell,2)
real(dp) :: cell_R(volume_mesh%ncell)

!Extract parameters
kd = 0.2d0
CFL = 0.9d0
conresbound = -12.0d0
niter = 15000

!Initialise
cell_phi(:) = 0.0d0 

!Set to initial guess 
! call guess_distance_field(cell_phi,BCactive,volume_mesh,meshdat)

!Solve
if (cm2dopt%dispt == 1) then
    write(*,'(A)')   '    +------------+------------+'
    write(*,'(A,A)') '    |    ittn    |','    dres    |'
end if
do solveiter=1,niter

    !Evaluate distance first gradient 
    phi_grad(:,:) = 0.0d0 
    do ee=1,volume_mesh%Nedge

        !Adjacent cells
        c3 = volume_mesh%edge(ee,3)
        c4 = volume_mesh%edge(ee,4)

        !Edge values 
        rc4 = cell_phi(c4)
        if (c3 .GT. 0) then !Cell value
            rc3 = cell_phi(c3)
        else !Boundary condition 
            if ((c3 .LT. 0) .AND. (BCactive(ee) == 1)) then    
                rc3 = 0.0d0 
            else
                rc3 = rc4
            end if 
        end if
        
        !Edge dimension properties
        edx = meshdat%edgedx(ee)
        edy = meshdat%edgedy(ee)

        !Edge distance (upwind)
        re = min(rc3,rc4)
        
        !Gradient flux
        gxf = edy*re
        gyf = -edx*re

        !Accumulate gradient
        phi_grad(c4,1) = phi_grad(c4,1) + gxf
        phi_grad(c4,2) = phi_grad(c4,2) + gyf
        if (c3 .GT. 0) then 
            phi_grad(c3,1) = phi_grad(c3,1) - gxf
            phi_grad(c3,2) = phi_grad(c3,2) - gyf
        end if 
    end do 
    phi_grad(:,1) = phi_grad(:,1)/meshdat%cell_area(:)
    phi_grad(:,2) = phi_grad(:,2)/meshdat%cell_area(:)
    phi_grad_mag(:) = (phi_grad(:,1)**2 + phi_grad(:,2)**2)

    !Evaluate distance second gradient
    phi_grad2(:) = 0.0d0 
    do ee=1,volume_mesh%Nedge

        !Adjacent cells
        c3 = volume_mesh%edge(ee,3)
        c4 = volume_mesh%edge(ee,4)
        
        !Edge values 
        rc4 = cell_phi(c4)
        gr4x = phi_grad(c4,1)
        gr4y = phi_grad(c4,2)
        if (c3 .GT. 0) then !Cell value
            rc3 = cell_phi(c3)
            gr3x = phi_grad(c3,1)
            gr3y = phi_grad(c3,2)
        else !Boundary condition 
            if ((c3 == -1) .AND. (BCactive(ee) == 1)) then 
                rc3 = 0.0d0 
            else
                rc3 = rc4
            end if 
            gr3x = gr4x
            gr3y = gr4y
        end if 

        !Vector between cell centres on this edge
        if (c3 .GT. 0) then 
            dfx = meshdat%cell_midp(c3,1) - meshdat%cell_midp(c4,1)
            dfy = meshdat%cell_midp(c3,2) - meshdat%cell_midp(c4,2)
        else
            dfx = meshdat%edge_midp(ee,1) - meshdat%cell_midp(c4,1)
            dfy = meshdat%edge_midp(ee,2) - meshdat%cell_midp(c4,2)
        end if 
        dfm2 = dfx**2 + dfy**2
        dfm = sqrt(dfm2)

        !Edge dimension properties
        edx = meshdat%edgedx(ee)
        edy = meshdat%edgedy(ee)
        
        !Edge FV gradient 
        gre_fvx = 0.5d0*(gr4x + gr3x)
        gre_fvy = 0.5d0*(gr4y + gr3y)

        !Edge fluxes
        flx_x = (((rc4 - rc3)/dfm)*(dfx/dfm) + gre_fvx)*0.5d0
        flx_y = (((rc4 - rc3)/dfm)*(dfy/dfm) + gre_fvy)*0.5d0 

        !Edge flux value 
        flx_edg = (flx_x*edy - flx_y*edx)

        !Accumulate to adjacent cells 
        phi_grad2(c4) = phi_grad2(c4) + flx_edg
        if (c3 .GT. 0) then 
            phi_grad2(c3) = phi_grad2(c3) - flx_edg
        end if 
    end do 
    phi_grad2(:) = phi_grad2(:)/meshdat%cell_area(:)

    !Residual 
    cell_R(:) = phi_grad_mag(:) - kd*phi_grad2(:) - 1.0d0 

    !Step   
    cell_phi(:) = cell_phi(:) - 0.5d0*CFL*meshdat%cell_area(:)*cell_R(:)

    !Display 
    dres = log10(sum(abs(cell_R(:)))/volume_mesh%ncell)
    ! dres = log10(sum(abs(cell_R(:))))
    !if (mod(solveiter,100) == 0) then 
        if (cm2dopt%dispt == 1) then
            write(*,'(A,I0,A,F9.5)') '        ',solveiter,'       ',dres
        end if 
    !end if 

    !Convergence check 
    if (dres .LE. conresbound) then 
        if (cm2dopt%dispt == 1) then
            write(*,'(A)') '           ** resconv **'
            write(*,'(A,I0,A,F9.5)') '        ',solveiter,'       ',dres
        end if 
        exit 
    end if 
end do 

!Store result
cell_d(:) = cell_phi(:)
grad_d(:,:) = -phi_grad(:,:)
grad2_d(:) = phi_grad2(:)
return 
end subroutine build_distance_field




!Subroutine to set distance field initial guess ================================================
subroutine guess_distance_field(cell_phi,BCactive,volume_mesh,meshdat)
implicit none 

!System data 
type(meshdata) :: meshdat 
type(vol_mesh_data) :: volume_mesh

!Variables - Import 
integer(in) :: BCactive(volume_mesh%nedge)
real(dp) :: cell_phi(volume_mesh%ncell)

!Variables - Local 
integer(in) :: ff,ee,cc,c3,c4,ctgt,cbase,nupdate,Nlevel
integer(in) :: cell_active(volume_mesh%ncell)
integer(in) :: cell_activeN(volume_mesh%ncell)
integer(in) :: cell_level(volume_mesh%ncell)
real(dp) :: sq_atotal,ncell,dist_current,dist_level

!Set zero distance where boundary conditions are active
cell_active(:) = 0
do ee=1,volume_mesh%Nedge
    c3 = volume_mesh%edge(ee,3)
    c4 = volume_mesh%edge(ee,4)
    if ((c3 == -1) .AND. (BCactive(ee) == 1)) then 
        cell_phi(c4) = 0.0d0 
        cell_active(c4) = 1
    end if 
end do
cell_activeN(:) = cell_active(:)

!Flood distance level from boundary conditions 
cell_level(:) = 0 
do ff=1,volume_mesh%Nedge

    !Add distance
    do ee=1,volume_mesh%Nedge

        !Adjacent cells
        c3 = volume_mesh%edge(ee,3)
        c4 = volume_mesh%edge(ee,4)

        !If active non active pair on internal edge
        if ((c3 .GT. 0) .AND. (c4 .GT. 0)) then 
            if (cell_active(c3) .NE. cell_active(c4)) then 
                if (cell_active(c3) == 0) then 
                    ctgt = c3 
                    cbase = c4 
                else
                    ctgt = c4 
                    cbase = c4 
                end if
                cell_level(ctgt) = ff 
                cell_activeN(ctgt) = 2
            end if 
        end if 
    end do 

    !Update active states
    nupdate = 0
    do cc=1,volume_mesh%ncell 
        if (cell_activeN(cc) == 2) then 
            cell_activeN(cc) = 1
            nupdate = nupdate + 1
        end if 
    end do 
    cell_active(:) = cell_activeN(:)

    !Exit with no cell updates
    if (nupdate == 0) then 
        exit 
    end if
end do 
Nlevel = maxval(cell_level(:))

!Set actual distance for each level 
dist_current = 0.0d0 
do ff=1,Nlevel

    !Find average length scale for this level 
    ncell = 0.0d0 
    sq_atotal = 0.0d0     
    do cc=1,volume_mesh%ncell 
        if (cell_level(cc) == ff) then 
            ncell = ncell + 1.0d0 
            sq_atotal = sq_atotal + sqrt(meshdat%cell_area(cc))
        end if 
    end do 

    !Set this levels distance 
    dist_level = (sq_atotal/ncell) + dist_current
    do cc=1,volume_mesh%ncell 
        if (cell_level(cc) == ff) then 
            cell_phi(cc) = dist_level
        end if 
    end do

    !Update current distance
    dist_current = dist_level
end do 
return 
end subroutine guess_distance_field




!Project to vertices function ================================================
function project_cell2vtx(cell_d,meshdat,volume_mesh,setsurfzero) result(vtx_d)
implicit none 

!Mesh data 
type(meshdata) :: meshdat 
type(vol_mesh_data) :: volume_mesh

!Variables - Import 
integer(in) :: setsurfzero
real(dp) :: vtx_d(volume_mesh%nvtx)
real(dp) :: cell_d(volume_mesh%ncell)

!Variables - Local 
integer(in) :: ii,cc,ee
real(dp) :: vtx_dN(volume_mesh%nvtx)
real(dp) :: dist,frac

!Assign vertex densities as weighted average of adjacent cells 
vtx_d(:) = 0.0d0 
do ii=1,volume_mesh%nvtx
    do cc=1,4
        if (meshdat%vtx_2_cell(ii,cc) .NE. 0) then 
            vtx_d(ii) = vtx_d(ii) + cell_d(meshdat%vtx_2_cell(ii,cc))*meshdat%cell_vtx_W(ii,cc)
        end if 
    end do     
end do 

!Force all surface vertices to zero distance 
if (setsurfzero == 1) then 
    do ee=1,volume_mesh%nedge
        if (volume_mesh%edge(ee,3) == -1) then 
            vtx_d(volume_mesh%edge(ee,1:2)) = 0.0d0 
        end if
    end do 
end if 

!Smooth vertex distances
do ii=1,volume_mesh%nvtx
    if (meshdat%vtx_vlnc(ii) == 2) then !Set valence 2 vertices to the distance weighted average of their endpoints
        dist = norm2(volume_mesh%vertices(meshdat%vtx_v2v(ii,2),:) - volume_mesh%vertices(meshdat%vtx_v2v(ii,1),:))
        frac = norm2(volume_mesh%vertices(ii,:) - volume_mesh%vertices(meshdat%vtx_v2v(ii,1),:))/dist
        vtx_dN(ii) = frac*vtx_d(meshdat%vtx_v2v(ii,2)) + (1.0d0 - frac)*vtx_d(meshdat%vtx_v2v(ii,1))
    else !Set as its base value
        vtx_dN(ii) = vtx_d(ii)
    end if
end do
vtx_d(:) = vtx_dN(:)
return 
end function project_cell2vtx


end module cellmesh2d_distF_mod