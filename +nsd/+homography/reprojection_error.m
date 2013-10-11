function [err] = reprojection_error(H, p1, p2)
%--------------------------------------------------------------------
% 
% Author: Jeffrey Byrne (jbyrne@ssci.com)
%
%--------------------------------------------------------------------


%% Input check
if ((size(p1, 1) ~= 3) | (size(p2, 1) ~= 3))
  error('Points must be homogeneous');
end
if size(p1, 2) ~= size(p2, 2)
  error('Inconsistent point correspondence');
end


%% Reprojection error (image to scene)
% p21 = H*p1; % image features reprojected to scene plane
% p21 = [p21(1,:)./p21(3,:); p21(2,:)./p21(3,:)];  % dehomogeneize
% p2 = p2(1:2,:);  % scene plane features
% err = sum((p2 - p21).^2, 1);


%% Reprojection error (scene to image)
p2r = H\p2; % scene plane reprojected to image plane 
p2r = [p2r(1,:)./p2r(3,:); p2r(2,:)./p2r(3,:)];  % dehomogeneize
p1 = p1(1:2,:);  % dehomogeneize
err = sum((p1 - p2r).^2, 1);


