function [inf_norms, eoc_inf, L2_norms, eoc_L2, H1_norms, eoc_H1] = test_solver()
   
   f = @sin;
   
   c1 = -(cos(1)*exp(3) - exp(1)*sin(1))/(2*(exp(4) + 1));
   c2 = (exp(1)*(cos(1) + exp(2)*sin(1)))/(2*(exp(4) + 1));
   sol = @(x) c1*exp(x) + c2*exp(-x) + sin(x)/2;
   solp = @(x) c1*exp(x) - c2*exp(-x) + cos(x)/2;
   
   iterations   = 10;
   inf_norms    = zeros(iterations,1);
   eoc_inf      = zeros(iterations,1);
   L2_norms     = zeros(iterations,1);
   eoc_L2       = zeros(iterations,1);
   H1_norms     = zeros(iterations,1);
   eoc_H1       = zeros(iterations,1);
   
   for i = 1:iterations
       n = 100 + (100*(i-1));
       h = 2/(n-1);
       v_h = signorini_solver(n, h, basis_quadrature(f, n, h));
       
%           meshw = 0.001;       % step size
%           X = -1:meshw:1;    % domain
%           f = sol(-X);      % range
%           Y = diff(f)/meshw;   % first derivative
%           f2 = fe_function(v_h,h,X);
%           Y2 = fe_function_prime(v_h,h,X);
%           plot(X(:,1:length(Y)),Y,X,f,X,f2,X,Y2)
       
       inf_norms(i) = inf_norm(v_h,h,sol);
       L2_norms(i)  = L2_norm(h,n,@(x) fe_function(v_h,h,x),sol);
       H1_norms(i)  = H1_norm(h,n,@(x) fe_function(v_h,h,x),sol,@(x) fe_function_prime(v_h,h,x),solp);
       if (i > 1)
           eoc_inf(i) = EOC(inf_norms(i), inf_norms(i-1), h, (2/(n-100-1)));
           eoc_L2(i) = EOC(L2_norms(i), L2_norms(i-1), h, (2/(n-100-1)));
           eoc_H1(i) = EOC(H1_norms(i), H1_norms(i-1), h, (2/(n-100-1)));
       end
       %, x, 0.25.*fe_phi(0,h,x), x, 0.25.*fe_phi(1,h,x)
       
   end
end

function [EOC] = EOC(norm, norm_prev, h, h_prev)
    EOC = (log(norm) - log(norm_prev)) / (log(h) - log(h_prev));
end

function [norm] = inf_norm(v_h, h, sol)
    inf_mesh = -1:h/1000:1;
    inf_distances = abs(sol(-inf_mesh) - fe_function(v_h,h,inf_mesh));
    norm = max(inf_distances);
end

function [norm] = L2_norm(h,n,f,sol)
    norm = sqrt(L2_norm_squared(h,n,f,sol,'L2'));
end

function [norm] = L2_norm_squared(h,n,f,sol,str)
    norm = 0;
    L2_func = @(x)(sol(-x) - f(x)).^2;
    for i=0:n-1
        t_im1 = max(-1, -1 + ((i-1)*h));
        t_i   = -1 + (i*h);
        t_ip1 = min(1, -1 + ((i+1)*h));
        if (t_im1 ~= t_i)
            norm = norm + quadgk(L2_func, t_im1, t_i);
        end
        if (t_ip1 ~= t_i)
            norm = norm + quadgk(L2_func, t_i, t_ip1);
        end
    end

end

function [norm] = H1_norm(h,n,f,sol,fp,solp)
    norm = L2_norm_squared(h,n,f,sol,'L2');
    norm = norm + L2_norm_squared(h,n,fp,@(x)-solp(x),'L2p');
    norm = sqrt(norm);
end
