function [] = set_paths()
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne <jebyrne@cis.upenn.edu>
%
%--------------------------------------------------------------------------


%% Disable name conflicts during
warning('off','MATLAB:dispatcher:nameConflict');  


%% Support paths
addpath(pwd);
addpath(fullfile(pwd, 'data'));
fprintf('[%s]: Adding path %s\n', mfilename, fullfile(pwd, 'data'));


%% Unpack Dependencies
deps = {'vlbenchmarks-1.0c','lightspeed-2.3.1','vlfeat-0.9.16','matlabPyrTools-1.3.1','export_fig','brisk-0.1.1','HCI-stereo','middlebury-flow'}; % pyrtools last
for k=1:length(deps)
  % Unpack
  depdir = fullfile(pwd,'deps',deps{k});
  if ~exist(depdir,'dir')
    fprintf('[%s]: Unpacking %s\n', mfilename, depdir);
    zipfile = strcat(depdir,'.zip');
    unzip(zipfile, depdir);
  end

  % Paths
  fprintf('[%s]: Adding path %s\n', mfilename, depdir');
  addpath(genpath(depdir),'-begin');
end


%% Compiling
if ~(exist('corrDn') == 3)
  origdir = pwd;
  cd(fullfile(origdir,'deps','matlabPyrTools-1.3.1','matlabPyrTools-1.3.1'));
  eval('compile_mex');
  cd(origdir);
end
if exist('vl_sift','file') ~= 3
  fprintf('[%s]: Compiling deps/vlfeat-0.9.16\n', mfilename');
  thisdir=pwd; cd(fullfile(thisdir,'deps','vlfeat-0.9.16','toolbox')); 
  vl_compile(false);  % ONLY ON WINDOWS, NEED TO RUN MAKE ON LINUX
  cd(thisdir);
end  
if ~(exist('repmat') == 3)
  origdir = pwd;
  cd(fullfile(origdir,'deps','lightspeed-2.3.1','lightspeed'));
  eval('install_lightspeed');
  cd(origdir);
end
if ~(exist('+nsd/+seedoflife/nesting_distance_mex') == 3)
  fprintf('[%s]: compiling +nsd/+seedoflife/nesting_distance_mex.cpp...\n', mfilename);
  origdir = pwd;
  cd(fullfile(origdir,'+nsd','+seedoflife'))
  mex('nesting_distance_mex.cpp')
  cd(origdir);
end


%% Initialization
if ~(exist('vl_kmeans') == 3)
  vl_setup();
end


%% Restore
warning('on','MATLAB:dispatcher:nameConflict');  

