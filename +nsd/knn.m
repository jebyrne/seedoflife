function [yh] = knn(X,y,Z,k)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Odd neighbors?
if ~nsd.util.isodd(k)
  warning('set k to be odd to break ties');
end


%% k-nearest neighbors
D = sqdist(X,Z);  % exhaustive pairwise distance
[D,k_sort] = sort(D,1,'ascend');  % nearest neighbors
yh = median(y(k_sort(1:k,:)));  % median label of k-nearest


%% TESTING
% >> nsd.knn(rand(4,8), rand(1,8), rand(4,2), 3)
