function [R] = rotation_matrix(r)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne 
% $Id: rotation_matrix.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------


%% Affine transformation matrix
if ~exist('r','var'), r  = 0; end
R = [cos(r) -sin(r); sin(r) cos(r)];  % 2D rotation matrix


