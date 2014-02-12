function [S] = saliency(im)

[d,di] = nsd.descriptor(im, fr);

D = nsd.seedoflife.reshape(d.^2,di);
for k=1:n_levels
  fr = nsd.detector.dense(im, 2^k, 2^k); 
  w = ndsum(D(:,:,k,:), 1:3);
  s(k) = nsd.util.imblur(w, 2^k);  
  S = S + s;
end

