function [] = test_solver()
   PHI = {};
   n = 10;
   h = 2/n;
   for i = 1:n
       PHI =  [PHI; {@(x){fe_phi(i,h,x)}}];
   end
   
   f = @sin;
   
   signorini_solver(n, basis_quadrature(f, PHI, n, h))
   
end

