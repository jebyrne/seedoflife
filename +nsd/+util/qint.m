function [p,y,a] = qint(ym1,y0,yp1) 
%QINT - quadratic interpolation of three adjacent samples
%
% [p,y,a] = qint(ym1,y0,yp1) 
%
% returns the extremum location p, height y, and half-curvature a
% of a parabolic fit through three points. 
% Parabola is given by y(x) = a*(x-p)^2+b, 
% where y(-1)=ym1, y(0)=y0, y(1)=yp1. 
% 
% https://ccrma.stanford.edu/~jos/sasp/Matlab_Parabolic_Peak_Interpolation.html

p = (yp1 - ym1)/(2*(2*y0 - yp1 - ym1)); 
y = y0 - 0.25*(ym1-yp1)*p;
a = 0.5*(ym1 - 2*y0 + yp1);

if isnan(p)
  p = 0;
end

%% References 
% FROM deps/vlfeat-0.9.16/vl/sift.c (for consistency)
%   /* quadratic interpolation */
%   double di = - 0.5 * (hp - hm) / (hp + hm - 2 * h0) ;
%   double th = 2 * VL_PI * (i + di + 0.5) / nbins ;
%   angles [ nangles++ ] = th ;
%
% Symbolic linear solution for least squares quadratic fit from three observations 
%   y(x) = A * x^2 + B * x + C
%   p = ( -B/2A , C - B^2/4A);  
