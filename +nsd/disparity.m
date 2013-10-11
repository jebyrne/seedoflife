function [imdisp, ij_left, ij_leftinright, ij_right, w_asgn, info] = disparity(imgin_left, imgin_right, opt)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne 
%
%--------------------------------------------------------------------------

%% Options
if ~exist('opt','var') || isempty(opt)
  opt = nsd.opts().disparity;
end


%% Correspondence - Epipolar
fprintf('[nsd.%s]: epipolar correspondence\n', mfilename);
[ij_left, ij_leftinright, ij_right, w_asgn, info] = nsd.correspondence(imgin_left, imgin_right, opt);


%% Sparse Disparity
d = ij_left(:,2) - ij_leftinright(:,2);  % positive disparity
d(d<0) = 0;  % remove geometrically invalid disparity
k_disp = sub2ind(size(info.f_ref.imgrey), round(ij_left(:,1)), round(ij_left(:,2))); 
imdisp = zeros(size(info.f_ref.imgrey));
imdisp(k_disp) = d(:);

