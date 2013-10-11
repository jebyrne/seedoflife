function [ij] = xy2ij(xy,d)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne 
% $Id: xy2ij.m 89 2012-09-16 23:41:37Z jebyrne $
%
%--------------------------------------------------------------------------

if ~exist('d','var')
  [xx,d] = min(size(xy));
end

ij = nsd.util.ij2xy(xy,d);  