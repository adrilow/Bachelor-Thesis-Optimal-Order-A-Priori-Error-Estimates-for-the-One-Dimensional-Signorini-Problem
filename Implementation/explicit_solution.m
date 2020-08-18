syms x c1 c2
u = symfun(c1*exp(x) + c2*exp(-x) + sin(x)/2, [x]);
up = symfun(c1*exp(x) - c2*exp(-x) + cos(x)/2, [x]);
eqns = [u(-1) >= 0, u(1) >=0, up(1)>=0, up(-1)<=0, u(1)*up(1)==0, u(-1)*up(-1)==0];
[sol1, sol2] = solve(eqns, [c1,c2], 'Real', true)

% sol1 = -(cos(1)*exp(3) - exp(1)*sin(1))/(2*(exp(4) + 1))
% sol2 = (exp(1)*(cos(1) + exp(2)*sin(1)))/(2*(exp(4) + 1))