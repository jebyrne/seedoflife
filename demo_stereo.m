function [] = demo_stereo(imleft, imright)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne 
%
%--------------------------------------------------------------------------
close all; 

%% Inputs
if nargin==0
  imleft = 'middlebury_stereo_teddy_im2.png';  % ./data
  imright = 'middlebury_stereo_teddy_im6.png';  % ./data
end


%% Options
opt = nsd.opts().disparity;
opt.do_debug = true;


%% Stereo!
imdisp = nsd.disparity(imleft,imright,opt);
nsd.show.disparity(imdisp);
