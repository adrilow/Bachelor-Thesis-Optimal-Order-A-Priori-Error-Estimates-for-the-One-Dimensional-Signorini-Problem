function [inf_norms, eoc_inf] = test_solver()
   
   
   f = @sin;
   
   c1 = -(cos(1)*exp(3) - exp(1)*sin(1))/(2*(exp(4) + 1));
   c2 = (exp(1)*(cos(1) + exp(2)*sin(1)))/(2*(exp(4) + 1));
   sol = @(x) c1*exp(x) + c2*exp(-x) + sin(x)/2;
   
   
   iterations = 10;
   inf_norms = zeros(iterations,1);
   eoc_inf = zeros(iterations,1);
   for i = 1:iterations
       n = 100 + (100*(i-1));
       h = 2/(n-1);
       v_h = signorini_solver(n, h, basis_quadrature(f, n, h));
       
       inf_mesh = -1:h/1000:1;
       inf_distances = abs(sol(-inf_mesh) - fe_function(v_h,h,inf_mesh));
       inf_norms(i) = max(inf_distances);
       
       if (i > 1)
           eoc_inf(i) = EOC(inf_norms(i), inf_norms(i-1), h, (2/(n-100-1)));
       end
       %, x, 0.25.*fe_phi(0,h,x), x, 0.25.*fe_phi(1,h,x)
       %plot(x, sol(-x), x, fe_function(v_h,h,x))
   end
end

function [EOC] = EOC(norm, norm_prev, h, h_prev)
    EOC = (log(norm) - log(norm_prev)) / (log(h) - log(h_prev));
end

