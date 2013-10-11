function [img] = imdisc(M,N,x,y,r,m,s)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
% $Id: imdisc.m 104 2013-02-04 16:45:20Z jebyrne $
%
%--------------------------------------------------------------------------

%% Inputs
if nargin == 0
  M = 128; N = 128;
  x = 64; y = 64;
  r = 32;
  m = 1;
  s = 1;
end
if isempty(x)
  x = round(N/2);
end
if isempty(y)
  y = round(M/2);
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
      if (floor(sqrt((i-x)^2+(j-y)^2)) <= r)
         mask(i,j) = 1;
      end
   end
end




