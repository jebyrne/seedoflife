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
deps = {'vlbenchmarks-1.0d','lightspeed-2.3.1','vlfeat-0.9.16','matlabPyrTools-1.3.1','export_fig-0.2','brisk-0.1.1','HCI-stereo','middlebury-flow'}; % pyrtools last
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

%% Compile MEX
try
  run('compile_mex.m');
catch ME
  fprintf('[%s]: compiling mex files failed - try running compile_mex.m to track down the problem', mfilename);
end


%% Restore
warning('on','MATLAB:dispatcher:nameConflict');  

