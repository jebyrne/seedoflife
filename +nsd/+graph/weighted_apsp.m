function [D,P] = weighted_apsp(G,W)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Inputs
[N,M] = graphmatch.get_size(G);
if (min(W) < 0)
  error('negative weights');
elseif length(W) ~= M
  error('invalid weight length');
end


%% Embedding distance (G)
[i,j,xx] = find(tril(G.adj,-1));  % zeros on main diagonal
E_sp = W(:) + 1E-6;  % minimum weight (base-M exponential encoding)
W_sp = sparse(i,j,E_sp,N,N); % edge weight matrix (column order of edge weights)
W_sp = max(W_sp,W_sp');  % symmetric
[i,j,W_sp] = find(W_sp);  % symmetric edge weights (column order)
[As,A,eil] = indexed_sparse(i,j,W_sp,N,N,struct('undirected',1));  % matlab-bgl


%% All pairs shortest path 
sp_opt.edge_weight = W_sp(eil);  % properly ordered edge weights
if nargout == 1
  sp_opt.algname = 'johnson';  % sparse, fast
  [D] = all_shortest_paths(As,sp_opt);  % matlab-bgl
else
  sp_opt.algname = 'floyd_warshall';  % slow, debugging paths
  [D,P] = all_shortest_paths(As,sp_opt);  % matlab-bgl
end


%% Debugging
% graphmatch.show_graph(G);
% p = graphmatch.path_from_apsp_pred(u,v,D,P);
% hold on; plot(G.xy(p,1),G.xy(p,2),'g.');



%% ARCHIVE
%{

%% Discretization: two levels only
nq = 100;
[N,N] = size(G.adj);
M = N*nq; % base-M
W_max = max(W); 
W_min = min(W);
W_scale = round((nq-1)*((W - W_min)./(W_max-W_min)))+1; % [W_min,W_max]->[1,nq]
Q(W_scale<(nq/2))  = W_scale(W_scale<(nq/2));  % 
Q(W_scale>=(nq/2)) = round(M*W_scale(W_scale>=(nq/2)));  % truncate

%% Decoding  
% k = find(D >= M);
% D(k) = fix(D(k)./M) + M*((D(k)./M)-(fix(D(k)./M)));  % decode
% D = ((W_max-W_min)*((D-1)./nq)) + W_min; % [1,nq]->[W_min,W_max]

%}