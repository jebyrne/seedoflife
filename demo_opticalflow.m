function [] = demo_opticalflow(im1, im2)
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

%% Examples
if nargin==0
  im1 = 'middlebury_flow_grove2_frame10.png';  % ./data
  im2 = 'middlebury_flow_grove2_frame11.png';  % ./data
end


%% Optical flow correspondence
[ij_1, ij_1in2, ij_2, w_asgn, info] = nsd.correspondence(im1, im2, nsd.opts().flow);  
imagesc(nsd.show.opticalflow(info.f_ref.imgrey, ij_1, ij_1in2));