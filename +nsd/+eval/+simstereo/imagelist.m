function [imglist] = imagelist(indir)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne 
%
%--------------------------------------------------------------------------


%% Inputs
%insubdir = fullfile(indir,'middlebury-stereo');
insubdir = indir;  % directory containing middlebury-stereo dataset


%% Image list
imglist(1).left = fullfile(insubdir, 'teddy', 'im2.png');
imglist(1).right = fullfile(insubdir, 'teddy', 'im6.png');
imglist(1).disp = fullfile(insubdir, 'teddy', 'disp2.png');  % disparity 2-6, relative to 2
imglist(1).dispscale = 4;
imglist(1).disc = fullfile(insubdir, 'teddy', 'disc2.png');  % discontinuitites
imglist(1).name = 'teddy';

imglist(2).left = fullfile(insubdir, 'cones', 'im2.png');
imglist(2).right = fullfile(insubdir, 'cones', 'im6.png');
imglist(2).disp = fullfile(insubdir, 'cones', 'disp2.png');
imglist(2).dispscale = 4;
imglist(2).disc = fullfile(insubdir, 'cones', 'disc2.png');  % discontinuitites
imglist(2).name = 'cones';

imglist(3).left = fullfile(insubdir, 'tsukuba', 'scene1.row3.col3.ppm');
imglist(3).right = fullfile(insubdir, 'tsukuba', 'scene1.row3.col5.ppm');
imglist(3).disp = fullfile(insubdir, 'tsukuba', 'truedisp.row3.col3.pgm');
imglist(3).dispscale = 8;
imglist(3).disc = fullfile(insubdir, 'tsukuba', 'discocc.row3.col3.png');  % discontinuitites
imglist(3).name = 'tsukuba';

imglist(4).left = fullfile(insubdir, 'venus', 'im2.ppm');
imglist(4).right = fullfile(insubdir, 'venus', 'im6.ppm');
imglist(4).disp = fullfile(insubdir, 'venus', 'disp2.pgm');
imglist(4).dispscale = 8;
imglist(4).disc = fullfile(insubdir, 'venus', 'disc2.png');  % discontinuitites
imglist(4).name = 'venus';

imglist(5).left = fullfile(insubdir, 'map', 'view1.png');
imglist(5).right = fullfile(insubdir, 'map', 'view5.png');
imglist(5).disp = fullfile(insubdir, 'map', 'disp1.png');
imglist(5).dispscale = 8;
imglist(5).disc = fullfile(insubdir, 'map', 'disc1.png');  % discontinuitites
imglist(5).name = 'map';

imglist(6).left = fullfile(insubdir, 'sawtooth', 'im2.ppm');
imglist(6).right = fullfile(insubdir, 'sawtooth', 'im6.ppm');
imglist(6).disp = fullfile(insubdir, 'sawtooth', 'disp2.pgm');
imglist(6).dispscale = 8;
imglist(6).disc = fullfile(insubdir, 'sawtooth', 'disc2.png');  % discontinuitites
imglist(6).name = 'sawtooth';

