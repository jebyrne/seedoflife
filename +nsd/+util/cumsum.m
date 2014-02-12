function [C] = cumsum(X)
%--------------------------------------------------------------------------
%
% Copyright (c) 2012 Jeffrey Byrne <jebyrne@cis.upenn.edu>
% $Id: cumsum.m 104 2013-02-04 16:45:20Z jebyrne $
%
%--------------------------------------------------------------------------

%% cumululative maximum over columns
C = X;
for k=2:size(X,1)
  C(k,:) = C(k-1,:) + X(k,:);
end