function [imflow] = opticalflow(imref, ij_ref, ij_refinobs, dminmax)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------


% Flow display (see flowToColor.m)
scale = 1;
u = ij_ref(:,1) - ij_refinobs(:,1);
v = ij_ref(:,2) - ij_refinobs(:,2);
%   k_invalid = find(sqrt(sum([u v].^2,2)) > 20);  % HACK: FIXME
%   u(k_invalid) = 0; % HACK
%   v(k_invalid) = 0; % HACK
rad = sqrt(u(:).^2 + v(:).^2);
rad = sort(rad,'ascend');
maxrad = max(-1, rad(end-round(length(rad)/10)));  % quantile
colorflow = computeColor(-(1/scale)*v/(maxrad+eps), -(1/scale)*u/(maxrad+eps));

[M,N] = size(imref);
k_flow = sub2ind([M N],round(ij_ref(:,1)),round(ij_ref(:,2)));
imflow = zeros(M*N,3);  % undefined is black
imflow(k_flow,:) = double(squeeze(colorflow));
imflow = mat2gray(reshape(imflow,[M N 3]));
