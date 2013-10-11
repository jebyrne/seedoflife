classdef nsdVlBenchmarkAdapter < localFeatures.GenericLocalFeatureExtractor & ...
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

% Authors: Karel Lenc, JEBYRNE

% AUTORIGHTS
  properties (SetAccess=public, GetAccess=public)
    Opts
    NsdArguments
    VlSiftArguments
    nsdopts
  end

  methods
    function obj = nsdVlBenchmarkAdapter(varargin)
      % def. arguments
      obj.Name = 'NSD';
      varargin = obj.configureLogger(obj.Name,varargin);
      obj.Opts = obj.checkInstall(varargin);
      if (~isempty(obj.Opts) && strcmp(obj.Opts{1},'Opts')) 
        obj.nsdopts = obj.Opts{2};
      else
        obj.nsdopts = nsd.opts();  % default
      end
      obj.ExtractsDescriptors = true;
      if isempty(obj.VlSiftArguments)
        obj.VlSiftArguments = {};
      end
      if isfield(obj.nsdopts,'vlsift')
        obj.VlSiftArguments = obj.nsdopts.vlsift; % HACK for oxford buildings
      end
        
      obj.UseCache = false; % JEBYRNE
    end

    function [frames descriptors info] = extractFeatures(obj, imagePath)
      import helpers.*;
      [frames descriptors] = obj.loadFeatures(imagePath,nargout > 1);
      if numel(frames) > 0; return; end;
      
      img = imread(imagePath);      
%      if(size(img,3)>1), img = mat2gray(rgb2gray(img)); end  % NEEDED only for consistency with demo_descriptors
      if(size(img,3)>1), img = (rgb2gray(img)); end  % needed for consistency with vlfeatsift detector      
           
      img = single(img); % If not already in uint8, then convert
      startTime = tic;
      if nargout == 1
        obj.info('Computing frames of image %s.',getFileName(imagePath));
        [frames] = vl_sift(img,obj.VlSiftArguments{:});
      else
        obj.info('Computing frames and descriptors of image %s.',...
          getFileName(imagePath));
        [frames] = vl_sift(img,obj.VlSiftArguments{:});
        [frames, descriptors, info] = obj.extractDescriptors(imagePath, frames);
      end
      
      timeElapsed = toc(startTime);
      obj.debug('%d Frames from image %s computed in %gs',...
        size(frames,2),getFileName(imagePath),timeElapsed);
      obj.storeFeatures(imagePath, frames, descriptors);
    end

    function [frames, descriptors, info] = extractDescriptors(obj, imagePath, frames)
      % extractDescriptor Extract SIFT descriptors of disc frames
      %   [DFRAMES DESCRIPTORS] = obj.extractDescriptor(IMG_PATH,
      %   FRAMES) Extracts SIFT descriptors DESCRIPTPORS of disc
      %   frames FRAMES from image defined by IMG_PATH. For the
      %   descriptor extraction, scale-space is used. Ellipses are
      %   converted to discs using their scale. The orientation of an
      %   oriented ellipse is dropped.
      import localFeatures.helpers.*;
      obj.info('Computing descriptors.');
      startTime = tic;
                
      if size(frames,1) > 4
        % Convert frames to disks
        frames = [frames(1,:); frames(2,:); getFrameScale(frames)];
      end
      if size(frames,1) < 4
        % When no orientation, compute upright SIFT descriptors
        frames = [frames; zeros(1,size(frames,2))];
      end

      % Preprocessing
      [img,xx,s] = nsd.preprocess(imagePath,obj.nsdopts.descriptor.pp);
      fr = frames; fr(2:-1:1,:) = (s*fr(1:2,:));  % xy -> ij
      if (s~=1)
        error('nonunity scale from preprocessing');  % idiot settings check
      end
      
      % Padding 
      if (min(size(img)) < 512)  % for small images in oxford buildings 
        [d,k] = min(size(img));  % minimum dimension for cspyr n_scales = 6
        if k(1) == 1
          img = padarray(img,[ceil((512 - d(1))/2) 0],'both','symmetric');
          di = ceil((512 - d(1))/2); dj = 0;
        else
          img = padarray(img,[0 ceil((512 - d(2))/2)],'both','symmetric');
          dj = ceil((512 - d(2))/2); di = 0;
        end
        fr(1:2,:) = fr(1:2,:) + repmat([di dj]',1,size(fr,2));  % offset frame for boundary padding
      end

      % Descriptors!
      [descriptors,di,fr,f] = nsd.descriptor(img,obj.nsdopts.descriptor,fr);
      elapsedTime = toc(startTime);
      obj.debug('Descriptors computed in %gs',elapsedTime);

  
      % Descriptor info
      info.di = di; info.fr = fr; info.f = f; % descriptor info            
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
