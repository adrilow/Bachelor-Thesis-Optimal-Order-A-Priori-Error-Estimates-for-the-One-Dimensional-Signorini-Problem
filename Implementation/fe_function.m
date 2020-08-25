function [y] = fe_function(v_h, h, x)
%UNTITLED Summary of this function goes here

y = zeros(size(x));

for i=1:max(size(v_h))
    y = y + (v_h(i) .* fe_phi(i-1,h,x));
end
end

