function [D] = sqdist(X,Y,Q,S)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne 
% $Id: opts.m 172 2013-08-05 21:10:12Z jebyrne $
%
%--------------------------------------------------------------------------

%% Inputs 
n_max = 8192;


%% All pairs distance
if max([size(X,2) size(Y,2)]) >= n_max
  % Break into pieces with sparsity hint  
  chunksize = 4096;
  [xx,k] = max([size(X,2) size(Y,2)]); if k(1)==1, chunkdim=1;, else, chunkdim=2; end
  if chunkdim == 1
    n_chunks = floor(size(X,2)/chunksize);
    fprintf('[nsd.util.sqdist]: breaking sqdist into %d chunks\n', n_chunks)
    for k=1:chunksize:(n_chunks-1)*chunksize
      fprintf('[nsd.util.sqdist][%d/%d]: chunked sqdist\n', k, (n_chunks-1)*chunksize)
      D(k:k+chunksize,:) = sqdist(X(:,k:k+chunksize),Y);
    end
    D(k+1:size(X,2),:) = sqdist(X(:,k+1:end),Y);  % remainder
  else
    n_chunks = floor(size(Y,2)/chunksize);
    fprintf('[nsd.util.sqdist]: breaking sqdist into %d chunks\n', n_chunks)
    for k=1:chunksize:(n_chunks-1)*chunksize
      fprintf('[nsd.util.sqdist][%d/%d]: chunked sqdist\n', k, (n_chunks-1)*chunksize)
      D(:,k:k+chunksize) = sqdist(X,Y(:,k:k+chunksize));
    end
    D(:,k+1:size(Y,2)) = sqdist(X,Y(:,k+1:end));  % remainder
  end      
elseif exist('pdist2') == 2
   D = pdist2(X',Y').^2;  % squared euclidean (stats toolbox, faster)
elseif exist('Q','var') && ~isempty(Q)
  D = sqdist(X,Y,Q);  % squared mahalanobis (lightspeed)
else
  D = sqdist(X,Y);  % squared euclidean (lightspeed)
end  

