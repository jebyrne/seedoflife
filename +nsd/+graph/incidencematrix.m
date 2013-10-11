function [A] = incidencematrix(g)
%--------------------------------------------------------------------------
%
% Author:  Jeffrey Byrne - jebyrne@cis.upenn.edu
% Copyright (c) 2009-2010 Jeffrey Byrne
%
%--------------------------------------------------------------------------

[V,E] = nsd.graph.size(g);
[I] = nsd.graph.node_index_matrix(g);
S = [ones(size(I,1),1); -ones(size(I,1),1)];
A = sparse([1:E 1:E]',I(:),S,E,V);

% if nargout > 1
%   Au = sparse([1:E],I(1:E),1,E,V);  % first endpoint only
%   Av = sparse([1:E],I(E+1:end),1,E,V);  % second endpoint only
% end

