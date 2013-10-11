function [opt] = opts()
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne 
%
%--------------------------------------------------------------------------


%% Global options
opt.verbose = true;


%% Preprocessing
opt.pp.maxsize = 128;  
opt.pp.rot = 0;
opt.pp.greyscale = 1;



%% Features
opt.features.spyr.n_bands = 8;  
opt.features.spyr.n_scales = floor(log2(opt.pp.maxsize))-2; 
opt.features.spyr.do_signed_orientation = false; 
opt.features.spyr.oe.max = inf;  % TESTING     
opt.features.spyr.oe.min = 0; 
opt.features.spyr.boundary = 'reflect1';
opt.features.pooling.mode = 'sum';
opt.features.pooling.m_support = [1 3 7]; 
opt.features.do_legacy_pooling = false;  % for old seedoflife only



%% Interest point detector
opt.detector.mode = 'canny';  
opt.detector.n_subsample = 1;
opt.detector.nested.r_nms = 5; % non-maximum suppression
opt.detector.nested.oemin = 1E-2; % non-maximum suppression
opt.detector.nested.n_bits = 6; % non-maximum suppression
opt.detector.canny.threshold = []; % empty for autothreshold
opt.detector.harris.radius = 2;
opt.detector.harris.threshold = 0.0005;
opt.detector.harris.sigma = 1;
opt.detector.features = opt.features;


%% Seed Of Life
sol.do_scale_invariant = false;  
sol.do_rotation_invariant = false;  
sol.mu = [];  % observation centroid for geometric regularization
sol.normalization.mode = 'none';  % see nsd.descriptor.seedoflife
sol.normalization.global = false;  
sol.normalization.e_robust_norm = 1E-6;
sol.do_scale_pooling = false;
sol.n_bits = inf;  % bits of precision
sol.n_lobes = 8;
sol.k_support = 3;
sol.dr = 3;
sol.r_nms = 3;
sol.do_affine = false;  % affine invariant interest points
sol.do_logspiral = true;  % difference normalization
sol.binarize = true;  % see nsd.descriptor.seedoflife

%% Flower of life


%% Descriptors
opt.descriptor.mode = 'seedoflife';
opt.descriptor.pp = opt.pp;  % preprocessing
opt.descriptor.nd = sol;  % BACKWARDS COMPATIBLE
opt.descriptor.sol = sol;  
opt.descriptor.detector = opt.detector;
opt.descriptor.features = opt.features;
opt.descriptor.verbose = opt.verbose;
opt.descriptor.sigma = 9;


%% Distance
opt.distance.mode = 'euclidean';
opt.distance.nesting.k_inlier = 0.7;  
opt.distance.nesting.k_outlier = 0.15;  
opt.distance.nesting.t_outlier = 0; 
opt.distance.r_gating = 16;
opt.distance.verbose = opt.verbose;


%% Assignment
opt.assignment.mode = 'greedy';
opt.assignment.grf.n_iter = 1000;


%% Refinement
opt.refinement.mode = 'none';
opt.refinement.r_adjacency = 8;


%% Correspondence
% opt.correspondence.assignment.mode = 'greedy';   % {'greedy','buttonhook'}
% opt.correspondence.refinement.mode = 'none';      % {'grf','gbp','none'}
% opt.correspondence.refinement.gbp.n_iter = 1000;
opt.correspondence.verbose = opt.verbose;
opt.correspondence.descriptor = opt.descriptor;
opt.correspondence.distance = opt.distance;  
opt.correspondence.assignment = opt.assignment;  
opt.correspondence.refinement = opt.refinement;  
%opt.correspondence.do_geometric_prior = false;  % geometric regularization
%opt.correspondence.gamma_geometric_prior = [0.01 0.01]./(opt.pp.maxsize/128);  % geometric regularization 


%% Disparity
opt.disparity = opt.correspondence;
opt.disparity.distance.mode = 'L2-scanline';
%opt.disparity.distance.mode = 'floweroflife-scanline';
opt.disparity.do_geometric_prior = false;


%% Optical flow
opt.flow = opt.correspondence;
opt.flow.distance.mode = 'L2-flow';
opt.flow.distance.r_gating = 16;


