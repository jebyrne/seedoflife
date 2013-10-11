function [K] = rips_complex(A, k)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------
% www.cs.dartmouth.edu/~afra/papers/smi10/vietoris.pdf


%% Inputs
if nargin == 1
  k = 2;
end
if k > 2
  error('unsupported k-simplex > 2');
end


%% Vietoris-Rips complex (K) from neighborhood graph (A)
[i,j] = find(tril(A,-1));  
K{1} = [1:size(A,1)];   % 0-simplex (oriented
K{2} = sort([i j]',1);  % 1-simplex (oriented)
for i=2:k
  X = K{i};  % i-simplexes
  Y = [];    % {i+1}-simplexes
  for j=1:size(X,2)  % for each i-simplex 
    tau = X(:,j);
    N = lower_neighbors(A,tau(1));
    for u=tau(2:end)'
      v = lower_neighbors(A,u);  
      N = intersect(N,v);
    end
    if ~isempty(N)
      for v=N'
        Y = [Y sort([tau; v])];
      end
    end
  end
  K{i+1} = Y;
end


%% ----
function N = lower_neighbors(A,u)
N = find(A(:,u));
N = N(N < u);
