function [f_v] = basis_quadrature(f, PHI, n)
%This function returns the vector of all (f, \phi_i)_{L^2}, for a given set
%of functions PHI = \{\phi_i | i = 1:n \}
    f_v = zeros(n)
    for i=1:n
        phi_i = PHI{i};
        f_v(i) = quadgk(@(x){f(x) * phi_i(x)}, -1,1);
end

