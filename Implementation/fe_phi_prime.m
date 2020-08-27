function [y] = fe_phi_prime(i,h,x)
%This is the i-th finite element basis function for a mesh width of h

t_i = -1 + (i*h);
t_im1 = -1 + ((i-1)*h);
t_ip1 = -1 + ((i+1)*h);

y = zeros(size(x));

ind1 = (x >= t_im1) & (x < t_i);
y(ind1) = 1/h;

ind2 = (x>= t_i) & (x < t_ip1);
y(ind2) = -1/h; 


end


