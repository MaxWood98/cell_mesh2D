%RBF Test
clearvars 

%Points
vtx = [0 0;
       1.0 1;
       1.0 -2; 
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



        % R(ii,jj) = abs(log(norm(vtx(ii,1) - vtx(jj,1))));
        % 
        % if ii == jj
        %     R(ii,jj) = 0;
        % end

        % if ii ~= jj 
        %     R(ii,jj) = max(norm(vtx(ii,1) - vtx(jj,1)),0.1);
        % else
        %     R(ii,jj) = 0;
        % end
        if R(ii,jj) > 0 && R(ii,jj) < mind
            mind = R(ii,jj); 
        end
    end
    mindist(ii) = mind;
end

Rdist = R;

maxdist = max(R,[],'all');
Rs = 50.5*maxdist;


%Smooth interpolation on nearby points and evaluate RBF

dmin = 0.1;
Wdpen = 0.1;

for ii=1:Npnts
    for jj=1:Npnts
        
        if ii ~= jj %Smooth
            if R(ii,jj) <= dmin*maxdist
                ii
                jj
                R(ii,jj) = R(ii,jj) + Wdpen*(R(ii,jj) - dmin*maxdist)^2;
            end
        end

        R(ii,jj) = wendlandc2(R(ii,jj),Rs);
    end
end


% midadj = 0.05;
% maxadj = 1 - midadj;
% for ii=1:Npnts
%     for jj=1:Npnts
%         if ii ~= jj 
%             if R(ii,jj) > maxadj
%                 R(ii,jj) = 0.5*(maxadj + R(ii,jj));
%             end
%         end
%     end
% end


% %Relax interpolation at interior points 
% Rsmooth = 0.01;%1e-1;
% for ii=2:Npnts-1
%     % R(ii,ii) = R(ii,ii) - Rsmooth;%/(mindist(ii)*Rs);
% 
%     % Rrow = R(ii,:);
%     % Rrow = Rrow(Rrow>0);
%     % min(Rrow)
%     % R(ii,ii) = R(ii,ii) - Rsmooth/(min(Rrow));
% 
%     % R(ii,ii) = 0.9;
% end

% Dw = 0.1;
% R(2,3) = R(2,3)*(1 - Dw/Rs);
% R(3,2) = R(3,2)*(1 - Dw/Rs);
% 
% 
% 
% dminfrac = 0.01;









% R(2,2) = 0.99;
% R(3,3) = 0.99;


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
        W = 0.0;
    else
        W = ((1 - d)^4)*(4*d + 1);
    end
end 