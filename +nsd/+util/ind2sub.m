function [i,] = ind2sub(matsize,k)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011-2012 Jeffrey Byrne
% $Id: ind2sub.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

%% Initialize
img_empty = zeros(imgsize);
img = zeros(imgsize);


%% Accumulate 
n_cols = size(X,2);
m_patch = round(sqrt(size(X,1)));
[i,j] = ind2sub(imgsize,1:n_cols);
for k=1:n_cols
  roi = reshape(X(:,k),[m_patch m_patch]);
  img = img + nsd.util.patch(img_empty,roi,[i(k) j(k)]);
end


