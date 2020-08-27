function [inf_norms, eoc_inf, L2_norms, eoc_L2, H1_norms, eoc_H1] = test_solver_L2()
   
   f = @(x) abs(x).^(-0.499);
   h_fine = 2^(-14);
   n_fine = 2^(14+1) + 1;
   [v_h_fine,B1_fine,B2_fine] = signorini_solver(n_fine, h_fine, basis_quadrature(f, n_fine, h_fine));
   sol = @(x) fe_function(v_h_fine, h_fine, x);
   solp = @(x) fe_function_prime(v_h_fine, h_fine, x);
   
   iterations   = 10;
   inf_norms    = zeros(iterations,1);
   eoc_inf      = zeros(iterations,1);
   L2_norms     = zeros(iterations,1);
   eoc_L2       = zeros(iterations,1);
   H1_norms     = zeros(iterations,1);
   eoc_H1       = zeros(iterations,1);
   
   for i = 1:iterations
       h = 2^(-i);
       n = 2^(i+1) + 1;
       h_prev = 2^(-(i-1));
       meshw = 1/2048;
       
       [v_h,~,~] = signorini_solver(n, h, basis_quadrature(f, n, h));
       v_h_refined = refine_v(h,h_fine,v_h);

       X1 = -1:meshw:1;    % domain
       f1 = sol(X1);      % range
       Y1 = diff(f1)/meshw;   % first derivative
       f2 = fe_function(v_h,h,X1);
       Y2 = fe_function_prime(v_h,h,X1);
       plot(X1(:,1:length(Y1)),Y1,X1,f1,X1,f2,X1(:,1:length(Y1)),Y2(:,1:length(Y1)));
       
       inf_norms(i) = inf_norm(v_h_refined,v_h_fine);
       L2_norms(i)  = L2_norm(v_h_refined,v_h_fine,B1_fine);
       H1_norms(i)  = L2_norm(v_h_refined,v_h_fine,B1_fine + B2_fine);
       if (i > 1)
           eoc_inf(i) = EOC(inf_norms(i), inf_norms(i-1), h, h_prev);
           eoc_L2(i) = EOC(L2_norms(i), L2_norms(i-1), h, h_prev);
           eoc_H1(i) = EOC(H1_norms(i), H1_norms(i-1), h, h_prev);
       end
       
   end
end

function [v_h_fine] = refine_v(h,h_fine,v_h)
    x = -1:h:1;
    xq = -1:h_fine:1;
    v_h_fine = interp1(x,v_h,xq)';
end

function [norm] = inf_norm(v_h_refined, v_h_fine)
    norm = max(abs(v_h_refined - v_h_fine));
end

function [norm] = L2_norm(v_h_refined, v_h_fine, B_1)
    dif = v_h_fine - v_h_refined;
    norm = sqrt(dif' * B_1 *dif);
end

