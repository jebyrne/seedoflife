function [X_max, X_maxindex] = maxfilter(X,m)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
% $Id: maxfilter.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

%% Input check
if ~nsd.util.isodd(m)
  error('max filter size must be odd');
end

%% Input parameters
[M,N] = size(X);


%% Max filter
X = padarray(X,[(m-1)/2 (m-1)/2],-1);  % zero padding
x = im2col(X,[m m],'sliding'); % sliding window
[x,k] = max(x,[],1);
X_max = reshape(x,[M N]);


%% Max filter index
if nargout == 2
  X_idx = reshape([1:M*N],M,N);  % linear index
  X_idx = padarray(X_idx,[(m-1)/2 (m-1)/2],-1);
  x_idx = im2col(X_idx,[m,m],'sliding');
  X_maxindex = reshape(x_idx(sub2ind(size(x_idx),k,1:length(k))), [M N]);
end
