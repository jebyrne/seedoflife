function [D] = distance(d_ref, di_ref, fr_ref, d_obs, di_obs, fr_obs, mode, opt)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne 
%
%--------------------------------------------------------------------------

%% Options
if ~exist('opt','var') || isempty(opt)
  opt = nsd.opts().distance;
end


%% Error check
if isstruct(d_ref) && isempty(strfind(mode,'floweroflife'))
  error('floweroflife requires floweroflife distance');
end


%% Pairwise Descriptor Distance
switch mode
    
  case {'L1'}
    if opt.verbose, fprintf('[nsd.%s]: L1 distance\n', mfilename); end
    D = vl_alldist2(d_ref,d_obs,'L1');  % ref -> obs 
    
  case {'L2-scanline'}
    if opt.verbose, fprintf('[nsd.%s]: L2-scanline distance\n', mfilename); end
    I=[]; J=[]; S=[];  % sparse index aggregation
    for k=1:max(fr_ref(1,:))
      k_ref = find(fr_ref(1,:)==k);  % reference scanline
      k_obs_scanline = find(abs(fr_obs(1,:) - k) <= 4);  % observed scanline bundle
      w = exp(-sqdist(d_ref(:,k_ref),d_obs(:,k_obs_scanline))./(32))';
      [X,Y] = meshgrid(k_ref,k_obs_scanline);
      I=[I;X(:)]; J=[J;Y(:)]; S=[S;w(:)];  % sparse index    
    end
    W = sparse(I,J,S,size(d_ref,2),size(d_obs,2));
    D = -W;  % sparse distance
    
  case {'L2-flow'}
    if opt.verbose, fprintf('[nsd.%s]: L2-flow gating distance (%1.2f px radius)\n', mfilename, opt.r_gating); end
    I=[]; J=[]; S=[];  % sparse index aggregation
    for k=1:size(d_ref,2)     
      k_box = find(sqdist(fr_obs(1:2,:),fr_ref(1:2,k)) <= (opt.r_gating.^2));  % gating
      w = exp(-sqdist(d_ref(:,k),d_obs(:,k_box))/(32))';     
      I=[I;k*ones(length(k_box),1)]; J=[J;k_box]; S=[S;w(:)];      
    end
    W = sparse(I,J,S,size(d_ref,2),size(d_obs,2));
    D = -W;  % sparse distance
    
  case {'nesting','nested'}
    if opt.verbose, fprintf('[nsd.%s]: nesting distance\n', mfilename); end
    [D] = nsd.seedoflife.nesting_distance_mex(single(d_ref),single(d_obs),di_ref.n_bands,di_ref.n_lobes,di_ref.n_scales,opt.nesting.k_inlier,opt.nesting.k_outlier,opt.nesting.t_outlier);
    %D(D>=opt.nesting.t_truncate) = NaN;
    D(D<0) = NaN;
            
  otherwise
    error('undefined distance mode ''%s''', mode);
end

