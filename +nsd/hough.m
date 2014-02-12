function [ij_max, v_max, Hq] = hough(W, ij_obs, dij_ref, houghsize, nq, opt)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Nomenclature
% W: n_obs x n_ref weight matrix
% ij_obs: ij position of observation interest points
% dij_ref: position of centroid relative to each reference interest point
% houghsize: unquantized dimensions of hough space (size of observed image)
% nq: quantization of hough space
% 
% ij_max: local of local maxima in hough space (sorted)
% v_max: value of global maximum in hough space (sorted)
% Hq: quantized hough voting 


%% Options
if ~exist('opt','var')
  opt.do_border = true;  % normalize by total number of possible votes
  opt.mn_filt = [2*nq 2*nq];  % final quantized interpolation filter size
  opt.var_filt = nq;  % final interpolation filter variance
end


%% Inputs
if nargin ~= 5
  error('invalid input')
end
M = houghsize(1); N = houghsize(2);  % hough space dimensions
[n_obs, n_ref] = size(W); % interest point distance


%% Hough voting (x,y)
ij_vote = zeros(n_obs*n_ref,2);
for j=1:n_ref
  imin = (j-1)*n_obs+1;
  imax = (j-1)*n_obs+n_obs;
  ij_vote(imin:imax,:) = ij_obs - repmat(dij_ref(j,:),n_obs,1);  % center vote
end
w_vote = W(:);
ij_vote = floor(ij_vote) + 1;
[k_valid] = nsd.util.inrect(ij_vote,[1 1 N M]);

% Quantization
H = sparse(ij_vote(k_valid,1),ij_vote(k_valid,2),w_vote(k_valid),M,N);
sumfun = @(x) sum(x(:))*ones(size(x));
%Hq = blkproc(H,[nq nq],round([nq/2 nq/2]),sumfun);  % full size, quantized
Hq = blkproc(H,[nq nq],sumfun);  % full size, quantized


%% Border reweight 
if opt.do_border 
  [V,U] = meshgrid(1:nq:N,1:nq:M); uv=[U(:) V(:)];
  n_center = size(uv,1);
  W_center = zeros(ceil(size(Hq)/nq));
  for k=1:n_center
    ij = repmat(uv(k,:),n_ref,1) + dij_ref;
    W_center(k) = n_ref ./ length(nsd.util.inmat(size(Hq), ij(:,1), ij(:,2)));
  end
  W_center(~isfinite(W_center)) = 0;
  W_center = imresize(W_center,size(Hq),'nearest');
  Hq = Hq .* W_center;  % border compensation
end


%% Interpolation
Hq = imfilter(Hq,fspecial('gaussian', opt.mn_filt, opt.var_filt)); 


%% Local maxima
[k_max,ij_max,v_max] = nsd.util.localmax(Hq);
[v_max,k] = sort(v_max,'descend');
ij_max = ij_max(k,:);
return;
 

