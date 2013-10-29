function [] = demo_correspondence(im1, im2)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne 
%
%--------------------------------------------------------------------------

%% Paths
if ~(exist('ndsum') == 2)
  fprintf('[%s]: running set paths\n', mfilename);
  run('set_paths.m');
end

%% Inputs
if nargin==0
  im1 = 'middlebury_stereo_teddy_im2.png';  % ./data
  im2 = 'middlebury_stereo_teddy_im6.png';  % ./data
end

%% Correspondence
[ij_1, ij_1in2, ij_2, w_asgn, info] = nsd.correspondence(im1, im2);  
nsd.show.matching(info.f_obs.imgrey, info.f_ref.imgrey, ij_1, ij_1in2, ij_2);
