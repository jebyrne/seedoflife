function [imgout] = imsmooth(img, s)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
% $Id: imsmooth.m 97 2012-12-24 16:17:52Z jebyrne $
%
%--------------------------------------------------------------------------

%% Inputs
if nargin == 1
  s = 1;
end

%% Blur
fsize = ceil(s*3) * 2 + 1;  % default size
h = fspecial('gaussian',[fsize fsize],s);
imgout = imfilter(double(img),h);