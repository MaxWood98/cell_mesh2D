%Gradient test

%Distance 
d = 0.9;

%Vertices
v1 = [0,0];
v2 = [1,0];

%Stepsize
h = 0.00001;

%Dv1_x
v1p = v1;
v1p(1) = v1p(1) + h;
v1m = v1;
v1m(1) = v1m(1) - h;
vip = get_vi(v1p,v2,d);
vim = get_vi(v1m,v2,d);
dv1_x = (vip - vim)/(2*h)
1 - d

%Function to evaluate point
function [vi] = get_vi(v1,v2,d)
    vi = v1 + d*(v2 - v1);
end