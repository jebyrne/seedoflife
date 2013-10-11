function [Y] = im2col(X,m,ij)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
% $Id: im2col.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

%% Input check
if ~nsd.util.isodd(m)
  error('sliding window size must be odd');
end


%% Input parameters
[M,N] = size(X);
mh = (m-1)/2;  % halfsize


%% Sliding window
if ~exist('ij','var')
  % slower
  X = padarray(X,[mh mh],0);  % zero padding
  Y = im2col(X,[m m],'sliding'); % sliding window
else
  % faster
  X = padarray(X,[mh mh],0);  % zero padding
  ij = ij + mh;  % zero padding offset
  [dj,di] = meshgrid(-mh:mh, -mh:mh);
  dj = dj(:); di = di(:);  
  n_pts = size(ij,1);
  i_slidingwin = zeros(m*m,n_pts);
  j_slidingwin = zeros(m*m,n_pts);
  for i=1:n_pts
    i_slidingwin(:,i) = ij(i,1) + di;
    j_slidingwin(:,i) = ij(i,2) + dj;
  end  
  k_slidingwin = sub2ind(size(X),i_slidingwin,j_slidingwin); 
  Y = X(k_slidingwin);
end


