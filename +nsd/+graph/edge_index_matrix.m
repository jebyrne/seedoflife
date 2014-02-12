function [I] = edge_index_matrix(A)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Edge index matrix I(i,j) = edge index for e(i,j)
N = size(A,1);
[i,j] = find(tril(A,-1));  % lower triangular only
k = [1:length(i)]';  % edge indexes from columnwise ordering
ii = [i;j]; jj = [j;i]; kk=[k;k]; % repeated symmetric edges
I = sparse(ii,jj,kk,N,N);  % construct sparse matrix

