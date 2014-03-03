function [imsalient, spyr] = saliency(im)
%--------------------------------------------------------------------------
%
% Copyright (c) 2014 Jeffrey Byrne
%
%--------------------------------------------------------------------------

nq = 4;
imgrey = nsd.util.im2gray(im);
fr = nsd.detector.dense(imgrey, nq, nq);

opt = nsd.opts();
opt.descriptor.verbose = false;
opt.descriptor.features.spyr.n_scales = floor(log2(min(size(imgrey))))-2;
opt.descriptor.sol.binarize = 0;
opt.descriptor.sol.do_logspiral = 1;
[d,di,fr,f] = nsd.descriptor(imgrey, opt.descriptor, fr);
W = sqrt((nsd.seedoflife.reshape(d,di)).^2);
W = min(W,0.5);  % truncate

% Pyramid decomposition
spyr = f.spyr.sepspyr;

% Low pass 
h = [0 -1 0; -1 4 -1; 0 -1 0]; 
spyr.lowpass = padarray(sqrt(conv2(spyr.lowpass, h, 'valid').^2), [1 1], 0, 'both');

% Bandpass
w = squeeze(mean(W, 2));
n_orientations = size(w,1);
for i=1:spyr.n_levels
  for j=1:spyr.n_basis
    spyr.bands{i,j} = zeros(size(spyr.bands{i,j}));
  end
  for j=1:n_orientations
    kappa = [spyr.filters.steer.inphase(-(j-1)*2*pi/n_orientations) 0];
    for k=1:spyr.n_basis
      spyr.bands{i,k} = real(spyr.bands{i,k}) + (kappa(k)*(1/(spyr.n_levels-i+1))*(imresize(full(sparse(floor(fr(1,:)/nq)+1, floor(fr(2,:)/nq)+1, squeeze(w(j,i,:)), floor(size(imgrey,1)/nq)+1, floor(size(imgrey,2)/nq)+1)), size(spyr.bands{i,k}), 'bicubic')));
    end
  end
end

% Filters
spyr.filters.f = abs(spyr.filters.f);
spyr.filters.f_order = real(spyr.filters.f_order);
spyr.filters.boundary = 0;

% Reconstruction!
imsalient = sepspyr.reconstruct(spyr);



