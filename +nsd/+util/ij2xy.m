function [xy] = ij2xy(ij,d)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne 
% $Id: ij2xy.m 89 2012-09-16 23:41:37Z jebyrne $
%
%--------------------------------------------------------------------------

if ~exist('d','var')
  [xx,d] = min(size(ij));
end

if d == 2
  % ij = [i(:) j(:)];
  xy = fliplr(ij);
else
  % ij = [i(:)'; j(:)'];
  xy = flipud(ij);
end