function [img] = centercrop(img, cropsize)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
% $Id: centercrop.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

%% Crop me
%rect = [XMIN YMIN WIDTH HEIGHT];
mu = size(img)/2;  % centroid
rect = round([mu(2)-(cropsize(2)/2) mu(1)-(cropsize(1)/2) cropsize(2) cropsize(1)]);
img = imcrop(img,rect);


