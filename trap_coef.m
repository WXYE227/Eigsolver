function [poles,weights] = trap_coef(center,r,d)
%generate coefficent for d-discrete, center = a, radius = b, circle

thetas = pi/d:(2*pi/d):(2*pi-pi/d);
poles =  r * exp(1i*thetas) + center;
weights =  r * exp(1i*thetas)/d;
end

