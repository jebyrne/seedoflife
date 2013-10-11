function [ijs,v,ccorr] = templatematch(img, tmpl)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
% $Id: templatematch.m 91 2012-12-12 17:03:07Z jebyrne $
%
%--------------------------------------------------------------------------

%% Options (FIXME)
opt.do_log = false;
opt.do_ncc = false;
opt.log_sigma = 2;
%opt.scales = [1:-0.125:0.25]';
opt.scales = 1;


%% Filtering: laplacian of gaussian
if opt.do_log 
  fsize = ceil(opt.log_sigma*3) * 2 + 1;  % default size
  h = fspecial('log',fsize,opt.log_sigma);
  img = imfilter(img,h,'same');
end


%% Cross correlation
scales = opt.scales;
n_scales = length(scales);
for k=1:n_scales
  tmplscale = imresize(tmpl,scales(k));
  if opt.do_log 
    tmplscale = imfilter(tmplscale,h,'same');
  end
  if opt.do_ncc
    y = normxcorr2(tmplscale, img);
    UV = (size(tmplscale)/2);
    ccorr(:,:,k) = y(round(UV(1):end-UV(1)),round(UV(2):end-UV(2)));  % 'same'
  else
    % cross correlation
    ccorr(:,:,k) = conv2(img,tmplscale,'same')./sqrt(sum(sum(tmplscale.^2)));
  end
end


%% Local maximum 
[ccorr,k_bestscale] = max(ccorr,[],3);  % max over scale
[k_max,ij_max,v_max] = nsd.util.localmax(ccorr); % local maximum
[v_max,k] = sort(v_max,'descend');  % sorted voting order
s_max = scales(k_bestscale(k_max)); 
ijs = [ij_max(k,:) s_max(k)];  % sorted voting order
ijs = ijs(1,:); % best only
v = v_max(1); % best only


%% Debugging
%figure(10); imagesc(ccorr); 

