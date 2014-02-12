function [D,P] = weighted_shortest_path(G,W,u,v)
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
if nargin < 4
  v = [];
end

%% Embedding distance (G)
[i,j,xx] = find(tril(G.adj,-1));  % zeros on main diagonal
E_sp = W(:) + 1E-6;  % minimum weight (base-M exponential encoding)
W_sp = sparse(i,j,E_sp,N,N); % edge weight matrix (column order of edge weights)
W_sp = max(W_sp,W_sp');  % symmetric
[i,j,W_sp] = find(W_sp);  % symmetric edge weights (column order)
[As,A,eil] = indexed_sparse(i,j,W_sp,N,N,struct('undirected',1));  % matlab-bgl


%% Single source shortest path (u->v)
sp_opt.edge_weight = W_sp(eil);  % properly ordered edge weights
if ~isempty(v)
  sp_opt.target = v;  % target node to terminate 
end
sp_opt.algname = 'dijkstra';  
[D,P] = shortest_paths(As,u,sp_opt);  % matlab-bgl


%% Path
if ~isempty(v)
  D = D(v);
  P = path_from_pred(P,v);
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