function [EOC] = EOC(norm, norm_prev, h, h_prev)
    EOC = (log(norm) - log(norm_prev)) / (log(h) - log(h_prev));
end
