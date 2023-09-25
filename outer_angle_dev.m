%Outer angle test
clc
cla reset
hold on

v1 = [0.1 0];
v2 = [-0.1 0];
vp = [0 1.0];

%Check convexity
ve1 = vp - v1;
ve2 = vp - v2;
[cp] = cp2d(ve1,ve2); 

%Check angle
cos_ang = dot(ve1,ve2)/(norm(ve1)*norm(ve2));
ang = acos(cos_ang);
ang = ang*(180/pi);
if cp < 0 
    ang = 360 - ang;
end
ang

%Plot
plot([v1(1) vp(1)],[v1(2) vp(2)],'r')
plot([v2(1) vp(1)],[v2(2) vp(2)],'b')

% axis tight
axis equal
% axis square
hold off


% dp = dot(v1,v2);
% sinang = cp2d(v1,v2)/(norm(v1)*norm(v2));
% 
% ang = asin(sinang)
% if dp < 0 
% 
% end
% 
% angd = (180/pi)*ang
% 
% 
%2d cross product
function [cp] = cp2d(v1,v2) 
    cp = v1(1)*v2(2) - v1(2)*v2(1);
end