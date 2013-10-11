classdef briskVlBenchmarkAdapter < localFeatures.GenericLocalFeatureExtractor & ...
    helpers.GenericInstaller
% localFeatures.VlFeatSift VlFeat vl_sift wrapper
%   localFeatures.VlFeatSift('OptionName',OptionValue,...) Creates new
%   object which wraps around VLFeat covariant image frames detector.
%   All given options defined in the constructor are passed directly
%   to the vl_sift function when called.
%
%   The options to the constructor are the same as that for vl_sift
%   See help vl_sift to see those options and their default values.
%
%   See also: vl_sift

% Authors: Karel Lenc, Jeffrey Byrne


% AUTORIGHTS
  properties (SetAccess=public, GetAccess=public)
    Opts
    BriskArguments
    VlSiftArguments
    nsdopts
  end

  methods
    function obj = briskVlBenchmarkAdapter(varargin)
      % def. arguments
      obj.Name = 'BRISK';
      varargin = obj.configureLogger(obj.Name,varargin);
      obj.Opts = obj.checkInstall(varargin);
      obj.ExtractsDescriptors = true;
      if isempty(obj.VlSiftArguments)
        obj.VlSiftArguments = {};
      end
      try, brisk('terminate'); end;
      brisk('init','threshold',60,'octaves',4);
      obj.UseCache = false;  % JEBYRNE

    end

    function [frames descriptors] = extractFeatures(obj, imagePath)
      import helpers.*;
      [frames descriptors] = obj.loadFeatures(imagePath,nargout > 1);
      if numel(frames) > 0; return; end;
      img = imread(imagePath);
      if(size(img,3)>1), img = rgb2gray(img); end
      img = single(img); % If not already in uint8, then convert
      
      obj.info('Computing BRISK descriptors of image %s.',getFileName(imagePath));
      startTime = tic;
      brisk('loadImage',fullfile('.',imagePath));      
      % detect keypoints
      frames=brisk('detect');
      % get the descriptors
      [frames descriptors]=brisk('describe');
      frames = frames';
      descriptors = descriptors';
      timeElapsed = toc(startTime);
      
      obj.debug('%d Frames from image %s computed in %gs',...
        size(frames,2),getFileName(imagePath),timeElapsed);
      obj.storeFeatures(imagePath, frames, descriptors);
    end

    function [frames descriptors] = extractDescriptors(obj, imagePath, frames)
      import localFeatures.helpers.*;
      obj.info('Computing descriptors.');
      startTime = tic;
          
      startTime = tic;      
      brisk('loadImage',fullfile('.',imagePath));      
      % detect keypoints
      frames=brisk('detect');
      % get the descriptors
      [frames descriptors]=brisk('describe');
      frames = frames';
      descriptors = descriptors';

      elapsedTime = toc(startTime);
      obj.debug('Descriptors computed in %gs',elapsedTime);
    end

    function sign = getSignature(obj)
      sign = [helpers.cell2str({rand()}),...
              helpers.cell2str(obj.VlSiftArguments)];
    end
  end

  methods (Access=protected)
    function deps = getDependencies(obj)
      deps = {helpers.VlFeatInstaller('0.9.14')};
    end
  end
end
