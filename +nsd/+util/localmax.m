function [k_max, ij_max, v_max, img] = localmax(A)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
% $Id: localmax.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

img = double(imregionalmax(A,8));  % localmax visualization
CC = bwconncomp(img);  % connected components
ij_max = zeros(CC.NumObjects,2);
for k=1:CC.NumObjects
  [i,j] = ind2sub(size(A),CC.PixelIdxList{k}');
  ij_max(k,:) = round(mean([i' j'],1));  % centroid
end
k_max = sub2ind(size(img),ij_max(:,1),ij_max(:,2)); % localmax index
v_max = A(k_max);  % value

% need nms suppression still single 8 connected is noisy

% Sorted
[v_max,k_sort] = sort(v_max,'descend');
k_max = k_max(k_sort);
ij_max = ij_max(k_sort,:);


