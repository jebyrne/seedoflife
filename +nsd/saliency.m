function [imsalient] = saliency(im)
%--------------------------------------------------------------------------
%
% Copyright (c) 2014 Jeffrey Byrne
%
%--------------------------------------------------------------------------

imgrey = nsd.util.im2gray(im);
fr = nsd.detector.dense(imgrey, 4, 4);

opt = nsd.opts();
opt.descriptor.features.spyr.n_scales = ceil(log2(min(size(imgrey))))-2;
%opt.descriptor.features.spyr.oe.max = 0.35;  % saturation
opt.descriptor.sol.binarize = 0;
opt.descriptor.sol.do_logspiral = 1;
[d,di] = nsd.descriptor(imgrey, opt.descriptor, fr);
W = (nsd.seedoflife.reshape(d,di).^2);

imsalient = zeros(size(imgrey));
for k=di.n_scales:-1:1
  s = ndsum(W(:,:,k,:), 1:3);  
  imsalient = imsalient + (2*(2^k + 1)*nsd.util.imblur(full(sparse(fr(1,:), fr(2,:), s, size(imgrey,1), size(imgrey,2))), 2*(2^k + 1)));
end




