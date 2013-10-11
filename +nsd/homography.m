function [H_lmse, lmse, ij_ref, ij_refinobs, ij_obs, f_ref, f_obs] = homography(imref, imobs, n_iter, do_debug)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne 
%
%--------------------------------------------------------------------------
% bagged least median of squares:
% http://research.microsoft.com/en-us/um/people/zhang/inria/publis/tutorial-estim/node25.html


%% Input check
if ~exist('n_iter','var') || isempty(n_iter)
  n_iter = 512; 
end
if ~exist('do_debug','var') || isempty(do_debug)
  do_debug = false; 
end


%% Nested descriptors
opt = nsd.opts().descriptor;
opt.do_geometric = false;
[d_ref,di_ref,fr_ref,f_ref] = nsd.descriptor(imref,opt);
[d_obs,di_obs,fr_obs,f_obs] = nsd.descriptor(imobs,opt);
[d_asgn,k_asgn] = min(sqdist(d_obs,d_ref),[],1);  % greedy assignment

ij_ref = fr_ref(1:2,:);
ij_obs = fr_obs(1:2,k_asgn);


%% Preconditioning
ijh_obs = nsd.util.homogenize(ij_obs);
K_obs = [1 0 0; 0 1 0; -mean(ij_obs,2)' 1]';  % translation
s = mean(sqrt(sum((nsd.util.dehomogenize(K_obs*ijh_obs)).^2,1)));  % mean deviation 
K_obs = [sqrt(2)/s 0 0; 0 sqrt(2)/s 0; -((sqrt(2)/s)*mean(ij_obs,2))' 1]';  % conditioning
ijhn_obs = K_obs*ijh_obs;  % normalization

ijh_ref = nsd.util.homogenize(ij_ref);
K_ref = [1 0 0; 0 1 0; -mean(ij_ref,2)' 1]';  % translation
s = mean(sqrt(sum((nsd.util.dehomogenize(K_ref*ijh_ref)).^2,1)));  % mean deviation 
K_ref = [sqrt(2)/s 0 0; 0 sqrt(2)/s 0; -((sqrt(2)/s)*mean(ij_ref,2))' 1]';  % conditioning
ijhn_ref = K_ref*ijh_ref;  % normalization


%% Bagged sample
[xx,k_sort] = sort(d_asgn,'ascend'); 
k_quantile = k_sort(1:round(0.75*length(k_sort)));  % percentile


%% Random sample concensus - Least median of squares
lmse = inf;  H_lmse = eye(3); % initial conditions
for i=1:n_iter
  % Random bagged sample of reference 
  k_homography = k_quantile(randperm(length(k_quantile),4));
  
  % Homography estimate on random sample
  [H] = planar_homography(ijhn_ref(:,k_homography), ijhn_obs(:,k_homography));

  % Median sum of square error of nested descriptors on test set
  if ~isempty(H) && (cond(H) < 1E1)
    ij_obs_test = nsd.util.dehomogenize(K_obs\(H*ijhn_ref));
    k_valid = nsd.util.inmat(size(f_obs.imgrey), ij_obs_test(1,:), ij_obs_test(2,:));
    k_valid = k_valid(1:1:end); % faster
    d_refinobs = nsd.seedoflife.seedoflife(f_obs, nsd.detector.frame(round(ij_obs_test(:,k_valid))',0,1), opt.nd);  % *not* rotation invariant frame!!! (TALK TO JEFF FOR FIX) 
    mse = median(sum((d_ref(:,k_valid) - d_refinobs).^2, 1));

    % Best?
    if mse < lmse
      lmse = mse;  % least median square error
      H_lmse = K_obs\(H*K_ref);  % best homography (inverse conditioning) 
      ij_refinobs = nsd.util.dehomogenize(K_obs\(H*ijhn_ref));
      if do_debug
        nsd.show.matching(f_obs.imgrey, f_ref.imgrey, ij_ref(:,k_valid)', ij_refinobs(:,k_valid)', ij_obs, figure(1));
      end
    end    
  end
end
return;



%% Planar homography - Direct Linear Transform
function [H] = planar_homography(p1, p2)
% Input check
if ((size(p1, 1) ~= 3) | (size(p2, 1) ~= 3))
  error('[sscv_planar_homography]: Points must be homogeneous');
end
if size(p1, 2) ~= size(p2, 2)
  error('[sscv_planar_homography]: Inconsistent point correspondence');
end

% 1. Compute a first approximation of the homography matrix
N = size(p1,2);
Xmatrix = zeros(9, 3*N);
for n = 1:N
  ss_p2 = [0 -1 p2(2,n); ...
             1 0 -p2(1,n); ...
             -p2(2,n) p2(1,n) 0];  % skew_symmetric
  a_n = kron(p1(:,n), ss_p2);
  Xmatrix(:, (n-1)*3+1:n*3) = a_n;
end
[U,S,V] = svd(Xmatrix');
HL = reshape(V(:,9), 3, 3);

% 2. Normalize by second largest singular value 
[U,S,V] = svd(HL);
H = HL ./ S(2,2);

% Choose the correct sign of H such that y'Hx > 0
xh = p1(:,1);
yh = p2(:,1);
if (yh'*H*xh < 0)
  H = -H;
end

