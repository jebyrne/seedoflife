function [X] = cummax(X)
%--------------------------------------------------------------------------
%
% Copyright (c) 2012 Jeffrey Byrne <jebyrne@cis.upenn.edu>
% $Id: cummax.m 165 2013-08-02 21:39:17Z jebyrne $
%
%--------------------------------------------------------------------------

%% cumululative maximum over rows
for k=2:size(X,1)
  X(k,:) = max(X(k-1,:),X(k,:));
end