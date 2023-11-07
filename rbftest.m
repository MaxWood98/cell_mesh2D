%RBF Test
clearvars 

%Points
vtx = [0 0;
       1 1;
       1.02 -2; 
       3 1; 
       5 0];

%Build dependance
[Npnts,~] = size(vtx);
R = zeros(Npnts,Npnts);
rbfloc = zeros(Npnts,1);
mindist = zeros(Npnts,1);
for ii=1:Npnts
    mind = 100000000;
    for jj=1:Npnts
        R(ii,jj) = norm(vtx(ii,1) - vtx(jj,1));
        if R(ii,jj) > 0 && R(ii,jj) < mind
            mind = R(ii,jj); 
        end
    end
    mindist(ii) = mind;
end



Rs = 50*max(R,[],'all');


for ii=1:Npnts
    for jj=1:Npnts
        R(ii,jj) = wendlandc2(R(ii,jj),Rs);
    end
end

%Relax interpolation at interior points 
Rsmooth = 1e-7;
for ii=2:Npnts-1
    R(ii,ii) = R(ii,ii) + Rsmooth/(mindist(ii)*Rs);
end



%Interpolate
Nint = 500;
% Ri = inv(R);
gamma = R\vtx(:,2);
% gamma = Ri*vtx(:,2);
xint = linspace(min(vtx(:,1)),max(vtx(:,1)),Nint)';
yint = zeros(Nint,1);
for ii=1:Nint
    for jj=1:Npnts
        rbfloc(jj) = abs(vtx(jj,1) - xint(ii,1));
        rbfloc(jj) = wendlandc2(rbfloc(jj),Rs);
    end 
    yint(ii,1) = sum(rbfloc.*gamma);
end



%% Plot

cla reset
hold on


plot(vtx(:,1),vtx(:,2),'r.','markersize',20)
plot(xint(:,1),yint(:,1),'b')

hold off
axis equal
% axis tight 
grid on


%% Functions

function [W] = wendlandc2(d,Rs)
    d = d/Rs;
    if d >= 1
        W = 0;
    else
        W = ((1 - d)^4)*(4*d + 1);
    end
end 