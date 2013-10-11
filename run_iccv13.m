function [] = run_iccv13(outdir, do_study)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne 
% $Id: run_iccv13.m 173 2013-08-10 01:31:12Z jebyrne $
%
%--------------------------------------------------------------------------

%% Inputs
if ~exist('outdir','var') || isempty(outdir)
  outdir = nsd.util.datefile();
end
nsd.util.mkdir(outdir);
basedir = pwd;


%% Datasets
middlebury_indir = './data/middlebury-stereo';


%% ICCV'13 plots
if ~exist('do_study')  || isempty(do_study)
  do_study = [1 1 1];
end


%% VGG Affine - SOL comparison
if do_study(1)
  cd(basedir);
  solopts = nsd.opts();
  solopts.descriptor.mode = 'seedoflife';
  solopts.descriptor.detector.n_subsample = 1;
  solopts.descriptor.pp.maxsize = [];
  solopts.descriptor.features.spyr.n_scales = 7;  
  solopts.descriptor.features.pooling.mode = 'sum';
  solopts.descriptor.features.spyr.boundary = 'reflect1';     
  solopts.descriptor.features.spyr.n_bands = 8;
  solopts.descriptor.features.spyr.n_scales = 7;
  solopts.descriptor.features.spyr.do_signed_orientation = false;  
  solopts.descriptor.sol.n_lobes = 8;
  solopts.descriptor.sol.do_affine = false;
  solopts.descriptor.sol.binarize = false;
  solopts.descriptor.sol.do_logspiral = true;
  solopts.vlsift = {'PeakThresh', 4};  % cannot use PeakThres=10, otherwise bark has one features.  http://www.vlfeat.org/benchmarks/overview/repeatability.html
  
  bsolopts = solopts;
  bsolopts.descriptor.sol.binarize = true;
  
  bsolDetector = nsd.eval.nsdVlBenchmarkAdapter('Opts',bsolopts);
  bsolDetector.Name = 'BinarySeedOfLife';  
  solDetector = nsd.eval.nsdVlBenchmarkAdapter('Opts',solopts);
  solDetector.Name = 'SeedOfLife';
  siftDetector = localFeatures.VlFeatSift('PeakThresh',4);  
  siftDetector.Name = 'SIFT';
  briskDetector = nsd.eval.briskVlBenchmarkAdapter();
  briskDetector.Name = 'BRISK';
  
  detectors = {bsolDetector, solDetector, siftDetector, briskDetector};
  datalist = {'bark', 'bikes', 'boat', 'graf', 'leuven', 'trees', 'ubc', 'wall'};
  %datalist = datalist([4 8 2 5 7 6 3 1]);  % include bark
  datalist = datalist([4 8 2 5 7 6 3]);   % leave out bark
  eval_vgg_affine(fullfile(outdir,'study_1'), detectors, datalist); 
end



%% VGG Affine - distance function comparison
if do_study(2)
  cd(basedir);
  solopts = nsd.opts();
  solopts.descriptor.mode = 'seedoflife';
  solopts.descriptor.detector.n_subsample = 1;
  solopts.descriptor.pp.maxsize = [];
  solopts.descriptor.features.spyr.n_scales = 7;  
  solopts.descriptor.features.pooling.mode = 'sum';
  solopts.descriptor.features.spyr.boundary = 'reflect1';   
  solopts.descriptor.features.spyr.n_bands = 8;
  solopts.descriptor.features.spyr.n_scales = 7;
  solopts.descriptor.features.spyr.do_signed_orientation = false;  
  solopts.descriptor.sol.n_lobes = 8;
  solopts.descriptor.sol.do_affine = false;
  solopts.descriptor.sol.binarize = false;
  solopts.descriptor.sol.do_logspiral = true;
  solopts.vlsift = {'PeakThresh',4};  % cannot use PeakThres=10, otherwise bark has one features.  http://www.vlfeat.org/benchmarks/overview/repeatability.html
  
  esolDetector = nsd.eval.nsdVlBenchmarkAdapter('Opts',solopts);
  esolDetector.Name = 'SeedOfLife (Euclidean Distance)';
  solDetector = nsd.eval.nsdVlBenchmarkAdapter('Opts',solopts);
  solDetector.Name = 'SeedOfLife (Nesting Distance)';  % WARNING: name must not change  
  
  detectors = {solDetector, esolDetector};
  datalist = {'bark', 'bikes', 'boat', 'graf', 'leuven', 'trees', 'ubc', 'wall'};
  datalist = datalist([4 8 2 5 7 6 3]);  % leave out bark 
  eval_vgg_affine(fullfile(outdir,'study_2'), detectors, datalist); 
end


%% Trade study
if do_study(3)
  cd(basedir);
  run_sol_tradestudy(middlebury_indir, fullfile(outdir,'study_3'), [10 20 30 50]);
end

