function [D] = reshape(d,di)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Reshape descriptor
n_desc = size(d,2);
D = reshape(d,[di.n_bands di.n_lobes di.n_scales n_desc]); 
