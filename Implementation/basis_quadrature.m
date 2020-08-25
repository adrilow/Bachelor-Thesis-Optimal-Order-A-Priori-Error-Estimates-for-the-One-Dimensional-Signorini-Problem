function [f_v] = basis_quadrature(f, n, h)
%This function returns the vector of all (f, \phi_i)_{L^2}, for a given set
%of functions PHI = \{\phi_i | i = 1:n \}
    f_v = zeros(n,1);
    for i=0:n-1
        t_im1 = max(-1, -1 + ((i-1)*h));
        t_ip1 = min(1, -1 + ((i+1)*h));
        f_v(i+1) = quadgk(@(x)f(x) .* fe_phi(i,h,x), t_im1, t_ip1);
    end
end

