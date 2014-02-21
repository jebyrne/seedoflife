function [d,di,fr,f,dm] = descriptor(imgfile, opt, fr)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne
%
%--------------------------------------------------------------------------


%% Options
if ~exist('opt','var') || isempty(opt);
  opt = nsd.opts().descriptor; 
end
if opt.verbose
  switch lower(opt.mode)
    case {'floweroflife'}
      fprintf('[nsd.%s]: descriptor=''%s'', detector=''%s''\n',mfilename,opt.mode,opt.detector.mode);
    case {'clover','cocentric','seedoflife','asol'}    
      if opt.sol.normalization.global
        str_norm = sprintf('%s-global',opt.sol.normalization.mode);
      else
        str_norm = sprintf('%s',opt.sol.normalization.mode);
      end
      if exist('fr','var') && ~isempty(fr)
        str_detector = 'frame';
      else
        str_detector = opt.detector.mode;
      end
      fprintf('[nsd.%s]: descriptor=''%s'', n_bands=%d, n_scales=%d, n_lobes=%d, binary=%d, rotsign=%d, logspiral=%d, norm=''%s'', pooling=''%s'', detector=''%s''\n',...
        mfilename,opt.mode,opt.features.spyr.n_bands, opt.features.spyr.n_scales, opt.sol.n_lobes,...
        opt.sol.binarize, logical(opt.features.spyr.do_signed_orientation>0),...
        opt.sol.do_logspiral, str_norm, opt.features.pooling.mode, str_detector);
    case {'sift'}
      fprintf('[nsd.%s]: descriptor=''%s'', detector=''%s''\n',mfilename,opt.mode,opt.detector.mode);
    otherwise
      error('invalid descriptor mode ''%s''', opt.mode);
  end
end 


%% Preprocessing
[im,imorig,imscale] = nsd.preprocess(imgfile, opt.pp);


%% Features
if ceil(log2(max(size(im)))) < (opt.features.spyr.n_scales+2)
  fprintf('[nsd.%s]: padding image to meet minimum scale requirements\n',mfilename);
  dmax = 2.^(opt.features.spyr.n_scales+2);
  im = padarray(im,max(ceil([(dmax - size(im,1))/2 (dmax - size(im,2))/2]),0),'both');
  if exist('fr','var') && ~isempty(fr)
    error('FIXME: padded image without frame shift');
  end
end

switch lower(opt.mode)
  case {'floweroflife'}
    f = nsd.seedoflife.features(im, opt.features);    
  case {'seedoflife','clover','cocentric','nsd','asol'}
    opt.features.do_legacy_pooling = true;  % FORCE ME
    f = nsd.seedoflife.features(im, opt.features);
  otherwise
    f = []; % no precomputed features 
end


%% Interest point detection (optional)
if ~exist('fr','var') || isempty(fr)
  [fr] = nsd.detector(im, opt.detector, f);  
end


%% Descriptor
opt.mode = lower(opt.mode);
switch opt.mode
  case {'asol'}
    [d,di] = nsd.seedoflife.asol(f, fr, opt.sol);

  case {'seedoflife'}
    [d,di] = nsd.seedoflife.seedoflife(f, fr, opt.sol);
    
  case {'clover','seedoflife-3'}
    opt.sol.k_support = 3;
    opt.sol.dr = 3;    
    [d,di,dd] = nsd.seedoflife.seedoflife(f, fr, opt.sol);
  case {'cocentric','seedoflife-0'}
    opt.sol.k_support = 2;
    opt.sol.dr = 0;
    opt.sol.n_lobes = 1;
    [d,di] = nsd.seedoflife.seedoflife(f, fr, opt.sol);        
  case 'sift'    
    imsift = im2single(mat2gray(im)); % image
    [xx,k_sort] = sort(fr(3,:),'ascend');  % ascending scale for vl_sift!
    fr_sift = [fr(2,k_sort); fr(1,k_sort); fr(3,k_sort); fr(4,k_sort)]; % ij -> xy
    [xx, d_sift] = vl_sift(imsift, 'frames', fr_sift);
    d_sift(:,k_sort) = d_sift;  % return to fr order
    d = nsd.util.normalize(double(d_sift));
    di = [];
    %imagesc(im); hold on; vl_plotsiftdescriptor(d_sift(:,1:10:end),fr_sift(:,1:10:end));
    
  otherwise
    error('invalid descriptor mode ''%s''', opt.mode);
end


%% Geometric regularization
% if opt.do_geometric
%   d = [d; nsd.descriptor.geometric(fr,f.center,opt.gamma)]; % centroid relative
%   di.is_geometric = true;
% end


%% Output
f.imorig = imorig;
f.imscale = imscale;
f.imgrey = im;

