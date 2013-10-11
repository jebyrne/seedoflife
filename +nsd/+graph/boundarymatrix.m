function [B] = boundarymatrix(K, k)
%--------------------------------------------------------------------------
%
% Author:  Jeffrey Byrne - jebyrne@cis.upenn.edu
% Copyright (c) 2009-2010 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Inputs
if ~iscell(K) && isnumeric(K)
  K = graphmatch.rips_complex(K,2);   % rips complex from input graph
elseif ~iscell(K) && (isstruct(K) && isfield(K,'adj'))
  K = graphmatch.rips_complex(K.adj,2);   % rips complex from input graph
elseif ~iscell(K)
  error('invalid rips complex K')
elseif k > length(K) 
  error('invalid k-boundary dimension');
end
if nargin == 1
  k = 1;
end

%% 0-boundary 
if k == 0
  X = K{1};
  [N] = size(X,2);
  B = eye(N,N);
  return;
end

%% Boundary homomorphism from oriented k-simplexes
X = K{k};
Y = K{k+1};
N = size(X,2);
M = size(Y,2);
B = sparse(N,M);     % Boundary matrix
for j=1:M
  f = Y(:,j)'; 
  for jj=1:length(f)
    fhat = [f(1:jj-1) f(jj+1:end)];
    i = ismember(X',fhat,'rows');  % ugh...
    B(i,j) = (-1)^(jj-1);
  end
end




%% References
% Hatcher: Algebraic Topology
% http://portal.acm.org/citation.cfm?id=1806721&dl=ACM&coll=DL&CFID=8295664&CFTOKEN=25761675
% http://compgeom.cs.uiuc.edu/~jeffe/pubs/pdf/surflow.pdf
% http://www.cs.brown.edu/courses/csci1490/slides/CS149-TUM.pdf
% http://www.acsu.buffalo.edu/~nagi/courses/684/Unimodularity.pdf
% http://knol.google.com/k/tony-illig/total-unimodularity/1e72jidh87s0i/4#
% perception.inrialpes.fr/~Boyer/ReadingGroup/GraphLaplacian-tutorial.pdf
