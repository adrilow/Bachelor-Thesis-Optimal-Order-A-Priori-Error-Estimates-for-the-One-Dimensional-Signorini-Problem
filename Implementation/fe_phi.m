function [y] = fe_phi(i,h,x)
%This is the i-th finite element basis function for a mesh width of h
t_i = -1 + (i*h);
t_im1 = -1 + ((i-1)*h);
t_ip1 = -1 + ((i+1)*h);

if (x >= t_im1 && x <= t_i)
    y = (x-t_im1)/h;
elseif (x>= t_i && x <= t_ip1)
    y = (t_ip1 - x)/h;
else
    y = 0;
end

