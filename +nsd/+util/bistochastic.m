function [What] = bistochastic(W) 
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
% $Id: bistochastic.m 203 2013-09-05 21:58:51Z jebyrne $
%
%--------------------------------------------------------------------------

%% Inputs
k_max = 15;
epsilon = 1E-2;
N = size(W,1);
M = size(W,2);
rowsum = 1;
colsum = N/M;

%% Bistochastic normalization
What = W;
for k=1:k_max
  % Sinkhorn normalization
  %What = What ./ repmat(sum(What,1),N,1);
  What = nsd.util.column_stochastic(What); % column
  %What = What ./ repmat(sum(What,2),1,M);
  What = nsd.util.column_stochastic(What')';  % row
  
  % Convergence?
  if all(abs(sum(What,1) - colsum) < epsilon) && all(abs(sum(What,2) - rowsum) < epsilon)
    break;
  end    
end

%% Converged?
if k == k_max
  warning('convergence failure - max error = %f', max(max(abs(sum(What,1) - colsum)),max(abs(sum(What,2) - rowsum))))
end


