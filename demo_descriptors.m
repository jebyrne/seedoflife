function [] = demo_descriptors(imreffile, imobsfile, k_demo, nsdopt)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne 
%
%--------------------------------------------------------------------------
close all;


%% Input
if ~exist('nsdopt','var') 
  nsdopt = nsd.opts();
end
if ~exist('imreffile','var') || isempty(imreffile)
  imreffile = 'peppers.png';
  %imreffile = 'cameraman.tif';
end
if ~exist('k_demo','var') || isempty(k_demo)
  k_demo = 3;
end
do_demo = zeros(1,100);
do_demo(k_demo) = 1;
imref = nsd.preprocess(imreffile, nsdopt.pp);


%% Observation
if ~exist('imobsfile','var') || isempty(imobsfile)
  fprintf('[%s]: Generating random homography\n', mfilename);
  %tij = [randi([-32 32],1); randi([-32 32],1)]; % central only
  %tij = [randi([-64 64],1); randi([-64 64],1)]; % central only
  %[A,A_prms] = nsd.util.affine_transform(tij, randn*pi/64, randn*0.2+1, randn*0.2+1,0*randn,0*randn); % random similarity
  %[A,A_prms] = nsd.util.similarity_transform([0 0], 0, 0.85); % random similarity
  %[A,A_prms] = nsd.util.similarity_transform([0 0], pi/4, 1); % random similarity
  %fprintf('[nsd.%s]: random affine transform (tij=[%1.1f,%1.1f], r=%1.2f, sx=%1.2f, sy=%1.2f, kx=%1.1f, ky=%1.1f)\n', mfilename, A_prms);
  
  %A = nsd.homography.random(128,128,4);
  [A,A_prms] = nsd.util.similarity_transform([0 0], -pi/3, 1); % random similarity
  imobs = nsd.util.imtransform(imref,A);
else
  imobs = nsd.preprocess(imobsfile, nsdopt.pp);
end


%% Nested shape descriptors
if do_demo(1)
  fprintf('[%s]: Demo 1 - visualization\n', mfilename);
  [d_ref,di_ref,fr_ref,f_ref] = nsd.descriptor(imref, nsdopt.descriptor);
  %nsd.show.descriptor(f_ref.imgrey,di_ref);
  nsd.show.frame(f_ref.imgrey,fr_ref);
  nsd.show.features(f_ref);
end


%% Nested shape descriptor - greedy matching 
if do_demo(2)
  fprintf('[%s]: Demo 2 - greedy matching \n', mfilename);
  [d_ref,di_ref,fr_ref] = nsd.descriptor(imref, nsdopt.descriptor);
  [d_obs,di_obs,fr_obs] = nsd.descriptor(imobs, nsdopt.descriptor);    
  [d_asgn,k_asgn] = min(sqdist(d_ref,d_obs),[],2);
  nsd.show.matching(imobs, imref, fr_ref(1:2,:)', fr_obs(1:2,k_asgn)', fr_obs(1:2,:)');  
end


%% Nested shape descriptor - correspondence
if do_demo(3)
  fprintf('[%s]: Demo 3 - correspondence \n', mfilename);
  [ij_ref,ij_refinobs,ij_obs,w,info] = nsd.correspondence(imref,imobs, nsdopt.correspondence);
  nsd.show.matching(info.f_obs.imgrey, info.f_ref.imgrey, ij_ref, ij_refinobs, ij_obs);  
  %nsd.show.distance(info.f_obs.imgrey, info.f_ref.imgrey, ij_refinobs, ij_ref);  
end


%% Nested shape descriptor - visualize matching descriptors
if do_demo(4)
  fprintf('[%s]: Demo 4 - visualize matching descriptors \n', mfilename);

  % Descriptors!
  [d_ref,di_ref,fr_ref,f_ref] = nsd.descriptor(imref,nsdopt.descriptor);
  [d_obs,di_obs,fr_obs,f_obs] = nsd.descriptor(imobs,nsdopt.descriptor);  
  
  % Choose matching interest points
  k_ref = nsd.show.nearest_interestpoint(f_ref.imgrey, fr_ref);
  k_obs = nsd.show.nearest_interestpoint(f_obs.imgrey, fr_obs);

  % Visualize
  D_ref = nsd.descriptor.reshape(d_ref,di_ref);
  D_ref = reshape(D_ref(:,:,:,k_ref), [di_ref.n_bands*di_ref.n_lobes di_ref.n_scales]);
  D_obs = nsd.descriptor.reshape(d_obs,di_obs);
  D_obs = reshape(D_obs(:,:,:,k_obs), [di_obs.n_bands*di_obs.n_lobes di_obs.n_scales]);  
  figure(1); clf; imagesc(D_ref,[0 nsdopt.descriptor.sol.features.spyr.oe.max]); title('ref');  colorbar;
  figure(2); clf; imagesc(D_obs,[0 nsdopt.descriptor.sol.features.spyr.oe.max]); title('obs');  colorbar;
end


%% NSD - Rotation invariance
if do_demo(7)
  fprintf('[%s]: Demo 7 - rotation invariance \n', mfilename);
  R = 0:(2*pi/64):2*pi;
  nsdopt.correspondence.descriptor.detector.mode = 'canny-rot';  % with rotation
  for r=R
    [A] = nsd.util.similarity_transform([0 0], -r, 1); % rotation
    imobs = nsd.util.imtransform(imref,A);
    [ij_ref,ij_refinobs,ij_obs,w,info] = nsd.correspondence(imref,imobs,nsdopt.correspondence);
    nsd.show.matching(info.f_obs.imgrey, info.f_ref.imgrey, ij_ref, ij_refinobs, ij_obs, figure(1));      
  end
end


%% NSD - Scale
if do_demo(8)
  fprintf('[%s]: Demo 8 - scale matching \n', mfilename);
  S = 1+ [0:0.125:1];
  for s=S
    [A] = nsd.util.similarity_transform([0 0], 0, s); % scale
    imobs = nsd.util.imtransform(imref,A);    
    [ij_ref,ij_refinobs,ij_obs,w,info] = nsd.correspondence(imref,imobs,nsdopt.correspondence);    
    nsd.show.matching(info.f_obs.imgrey, info.f_ref.imgrey, ij_ref, ij_refinobs, ij_obs, figure(2));      
  end
end


%% Covariance
if do_demo(9)
  fprintf('[%s]: Demo 9 - covariance \n', mfilename);
  if nsdopt.descriptor.sol.do_logspiral_difference 
    warning('logspiral difference enabled - forcing disabled');
    nsdopt.descriptor.sol.do_logspiral = false;
  end
  nsdopt.descriptor.detector.n_subsample = 8;  % for visualization
  [d,di,fr,f] = nsd.descriptor(imref, nsdopt.descriptor);
  [E,Ec] = nsd.descriptor.covariance(d,di,fr); 
  nsd.show.covariance(f.imgrey,fr(1:2,:)',E);
end



%% Nested shape descriptor - assignment
if do_demo(10)
  fprintf('[%s]: Demo 10 - assignment \n', mfilename);
  [d_ref,di_ref,fr_ref,f_ref] = nsd.descriptor(imref, nsdopt.descriptor);
  [d_obs,di_obs,fr_obs,f_obs] = nsd.descriptor(imobs, nsdopt.descriptor);    

  A = nsd.descriptor.assignment(d_ref,di_ref,d_obs,di_obs);
  
  while(1)
    k_ref = nsd.show.nearest_interestpoint(f_ref.imgrey, fr_ref);    
    figure(10); subplot(1,2,1); imagesc(f_ref.imgrey); colormap(gray); axis image;
    hold on; plot(fr_ref(2,k_ref),fr_ref(1,k_ref),'g.'); hold off;
    figure(10); subplot(1,2,2); imagesc(f_obs.imgrey); colormap(gray); axis image;
    hold on; plot(fr_obs(2,A(k_ref,:)),fr_obs(1,A(k_ref,:)),'g.'); hold off;
    drawnow;
  end
%  nsd.show.matching(imobs, imref, fr_obs(1:2,:)', fr_ref(1:2,k_asgn)');  
end


%% Flower of life
if do_demo(11)
  fprintf('[%s]: Demo 11 - flower of life \n', mfilename);
  [d_ref,di_ref,fr_ref,f_ref] = nsd.descriptor(imref, nsdopt.descriptor);
  [d_obs,di_obs,fr_obs,f_obs] = nsd.descriptor(imobs, nsdopt.descriptor);        
%  [U] = nsd.floweroflife.descriptor(fr_ref,8,32);  
%  [V] = nsd.floweroflife.descriptor(fr_obs,8,32);  
  
  keyboard
%   for k=1:4
%     for j=1:4
%       k_sel = find(U{j,k}(:,k_ref));
%       nsd.show.frame(f_ref.imgrey, fr_ref(:,k_sel), figure(1)); hold on
%       plot(fr_ref(2,k_ref),fr_ref(1,k_ref),'r.','MarkerSize',10);
%       U{j,k}(k_sel,k_ref)
% 
%       k_sel = find(V{j,k}(:,k_obs));
%       nsd.show.frame(f_obs.imgrey, fr_obs(:,k_sel), figure(2)); hold on
%       plot(fr_obs(2,k_obs),fr_obs(1,k_obs),'r.','MarkerSize',10);
%       V{j,k}(k_sel,k_obs)
%       
%       pause;
%     end
%   end
%   keyboard
    %k_ref = nsd.show.nearest_interestpoint(f_obs.imgrey, fr_obs)

  while(1)
    k_ref = nsd.show.nearest_interestpoint(f_ref.imgrey, fr_ref);
    %k_ref = 89;
    for k=size(U,2):-1:1
      u = U(:,k:end);
      v = V(:,k:end);
      D = nsd.floweroflife.distance(u,v,d_ref,d_obs,[],[]);  
      figure(2); plot(full(-D(k_ref,:)));  
      w = full(-D(k_ref,:)); k_valid = find(isfinite(w));
      w = exp((w(k_valid) - mean(w(k_valid))));
      nsd.show.imscatter(imobs,fr_obs(1:2,k_valid)', w, figure(1), 100);
      pause
    end
  end
  
  while(1)
    k_ref = nsd.show.nearest_interestpoint(f_ref.imgrey, fr_ref);
    D = nsd.floweroflife.distance(U,V,d_ref,d_obs,[],[]);
    figure(2); plot(full(-D(k_ref,:)));  
    w = full(-D(k_ref,:)); k_valid = find(isfinite(w));
    w = exp((w(k_valid) - mean(w(k_valid))));
    nsd.show.imscatter(imobs,fr_obs(1:2,k_valid)', w, figure(1), 100);
    pause
  end
    return;

%   for k=1:7
%     k_sel = find(U{k}(:,k_ref));
%     clf; nsd.show.frame(f_ref.imgrey, fr_ref(:,k_sel), figure(3)); hold on
%     plot(fr_ref(2,k_ref),fr_ref(1,k_ref),'r.','MarkerSize',10);
%     pause;
%   end
%   return;
  
  
%  keyboard
  %[d_asgn,k_asgn] = min(D,[],2);
  [u,v] = find(D);
  [dists, perm] = sort(nonzeros(D(:)),'ascend');
  %[aIdx bIdx] = ind2sub([size(d_ref,2), size(d_obs,2)], perm(1:nnz(D)));
  k_asgn = benchmarks.helpers.greedyBipartiteMatching(size(d_ref,2), size(d_obs,2), [u(perm) v(perm)]);  % vlbenchmarks
   k_valid = find(k_asgn);
%   %[k_asgn] = nsd.floweroflife.assignment(U,-D,V);
   nsd.show.matching(imobs, imref, fr_ref(1:2,k_valid)', fr_obs(1:2,k_asgn(k_valid))', fr_obs(1:2,k_asgn(k_valid))');
end


%% Debugging
if do_demo(12)
  fprintf('[%s]: Demo 12 - NSD debugging\n', mfilename);
  [ij_ref, ij_refinobs, ij_obs, w_asgn, info] = nsd.correspondence(imref, imobs, nsdopt.correspondence);  
  nsd.show.distance(info.f_obs.imgrey, info.f_ref.imgrey, ij_ref, ij_refinobs, ij_obs, ...
    info.d_obs,info.d_ref,info.di_obs,info.di_ref,info.fr_obs,info.fr_ref, ...
    info.f_obs,info.f_ref,info.D);
end
