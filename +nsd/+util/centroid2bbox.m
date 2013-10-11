function [bbox] = centroid2bbox(ij,matsize)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
% $Id: centroid2bbox.m 90 2012-10-02 15:07:56Z jebyrne $
%
%--------------------------------------------------------------------------

% bb = [x y width height]
bbox(1) = ij(2)-(matsize(2)/2);
bbox(2) = ij(1)-(matsize(1)/2);
bbox(3) = matsize(2);
bbox(4) = matsize(1);
bbox = round(bbox);

