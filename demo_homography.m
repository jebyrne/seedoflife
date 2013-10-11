function [] = demo_homography(tmplfile)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne 
%
%--------------------------------------------------------------------------
close all;

%% Inputs
if nargin == 0
  tmplfile = 'peppers.png';
end


%% Preprocessing
imtmpl = nsd.preprocess(tmplfile);
imtmpl = padarray(imtmpl,[48 48],0,'both');
imtmpl = nsd.preprocess(imtmpl,nsd.opts().pp,true);


%% Random homography
fprintf('[nsd.%s]: random homography \n', mfilename);
[A] = nsd.homography.random(size(imtmpl,1),size(imtmpl,2),8);
imobs = nsd.util.imtransform(imtmpl,A);


%% Robust homography estimation
[H,mse,ij_tmpl,ij_tmplinobs,ij_obs,f_tmpl,f_obs] = nsd.homography(imtmpl, imobs, 256, true);
nsd.show.matching(f_obs.imgrey, f_tmpl.imgrey, ij_tmpl', ij_tmplinobs', ij_obs',[]); 
fprintf('[%s]: mse = %f\n', mfilename, mse);


