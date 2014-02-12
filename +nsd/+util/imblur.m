function [imgout,h] = imblur(img, s)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
% $Id: imblur.m 149 2013-03-16 02:34:13Z jebyrne $
%
%--------------------------------------------------------------------------

%% Inputs
if ~exist('s','var') || isempty(s)
  s = 1;
end

%% Blur
fsize = ceil(s*3) * 2 + 1;  % default size
h = fspecial('gaussian',[fsize fsize],s);
imgout = imfilter(double(img),h);