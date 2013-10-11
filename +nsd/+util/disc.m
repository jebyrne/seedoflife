function [img] = disc(M,N,x,y,r,m,s)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
% $Id: disc.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

error('replaced by imdisc');

%% Inputs
if nargin == 0
  M = 128; N = 128;
  x = 64; y = 64;
  r = 32;
  m = 1;
  s = 1;
end

%% Outputs
img = zeros(M,N);

%% Draw circle at x,y,r with magnitude m
for k=1:length(x)
  img_k = circMask(x(k),y(k),r(k),M,N);
  img(find(img_k)) = m(k);
end

%% Blur
if s > 0
  fsize = ceil(s*3) * 2 + 1;  % default size
  h = fspecial('gaussian',[fsize fsize],s);
  img = imfilter(img,h);
end

function [mask] = circMask(x,y,r,X,Y)
mask = zeros(X,Y);
for i=1:X
   for j=1:Y
      if (sqrt((i-x)^2+(j-y)^2) <= r)
         mask(i,j) = 1;
      end
   end
end




