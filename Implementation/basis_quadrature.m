function [f_v] = basis_quadrature(f, PHI, n, h)
%This function returns the vector of all (f, \phi_i)_{L^2}, for a given set
%of functions PHI = \{\phi_i | i = 1:n \}
    f_v = zeros(n);
    for i=1:n
        phi_i = PHI{i};
        t_im1 = max(-1, -1 + ((i-1)*h));
        t_ip1 = min(1, -1 + ((i+1)*h));
        f_v(i) = quadgk(@(x){f(x) * phi_i(x)}, t_im1, t_ip1);
    end
end

