function [h] = Adapted_h(sigma)
% sigma: noise standard deviation
% h: adapted smoothing parameter

if sigma<5 
    h=0.8*sigma; 
elseif sigma<15
    h=0.6*sigma;
elseif sigma<25
    h=0.3*sigma; 
elseif sigma<35
    h=0.3*sigma;
elseif sigma<45
    h=0.2*sigma; 
else h=0.15*sigma; 
end

end