function [k_nms] = nms(W, ij, r)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
% $Id: nms.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------


%% Preprocessing: consider W at ij only
k = sub2ind(size(W),ij(1,:),ij(2,:));
w = W(k);
W = zeros(size(W));
W(k) = w;


%% Max filter 
[W_max, W_maxindex] = nsd.util.maxfilter(W,r);


%% Non-maximum suppression
w = W(k);
w_max = W_max(k);
k_nms = find(w == w_max);

