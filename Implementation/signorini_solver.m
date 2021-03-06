function [v_h, B_1, B_2] = signorini_solver(n, h, f)
%This function solves the 1D-Signorini Problem on the interval (-1,1)
%by a Finite Elements approximation with a mesh width of h
%   The actual problem solved is the following:
% $$\min \{ F(v_h) \,|\, v_h \in V_h \, v_h(-1) \geq 0,\, v_h(1) \geq 0 \}
% F(v_h) = \int_{-1}^1\frac{1}{2} (v_h')^2 + \frac{1}{2} v_h^2 - fv_h dx$$
%   The solution is the coefficient vector of the standard basis
%   representation of the finite element solution function v_h.
%   The parameter n controls the mesh width
%   The function f is the force. A quadrature of appropriate order
%   will be used for the rightmost part of the functional.

    B1_diag = (2/3) * ones(n,1);
    B1_diag(1) = 1/3; B1_diag(n) = 1/3;
    B1_codiag = (1/6) * ones(n,1);

    B2_diag = 2 * ones(n,1);
    B2_diag(1) = 1; B2_diag(n) = 1;
    B2_codiag = -1 * ones(n,1);

    B_1 = h     * spdiags([B1_codiag B1_diag B1_codiag], -1:1, n, n);
    B_2 = (1/h) * spdiags([B2_codiag B2_diag B2_codiag], -1:1, n, n);


    A_diag = zeros(n);
    A_diag(1) = -1; A_diag(n) = -1;
    b = zeros(n,1);
    A = spdiags(A_diag, 0, n, n);

    options = optimoptions('quadprog', 'Display', 'off');
    v_h = quadprog(B_1+B_2, -f, A, b, [],[],[],[],[],options);
end
