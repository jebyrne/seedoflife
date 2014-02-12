function [B] = repshift(A,p)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
% $Id: repshift.m 120 2013-02-25 19:35:21Z jebyrne $
%
%--------------------------------------------------------------------------

%% Shift array with repeated padding (see circshift docs)
idx = cell(1, ndims(A));
for k = 1:ndims(A)
    m = size(A,k);
    i = (0:m-1)-p(k);
    idx{k} = min(max(i,0),(m-1)) + 1;
end
B = A(idx{:});

