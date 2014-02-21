function [] = demo_hough(tmplfile, bgfile)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne 
%
%--------------------------------------------------------------------------
close all;

%% Paths
if ~(exist('ndsum') == 2)
  fprintf('[%s]: running set paths\n', mfilename);
  run('set_paths.m');
end

%% Inputs
if nargin == 0
  tmplfile = fullfile('data','vla.png');
  bgfile = fullfile('data','nimitz.jpg');

  tmpl = nsd.preprocess(tmplfile);
  tmpl = padarray(tmpl, [32 32], 0);
  img_bg = imresize(nsd.preprocess(bgfile), size(tmpl));
  img = imresize(nsd.util.im2gray(bgfile), 0.5);

  % Random transformation
  tij = [randi([-32 32],1); randi([-32 32],1)]; % central only
  r = rand*(pi/4) - pi/8;
  s = randn*0.2+1;
  [A,A_prms] = nsd.util.affine_transform(tij,r,s,s,0,0); % random similarity
  img = nsd.util.imtransform(tmpl,A);
  img = 0.5*img + 0.5*img_bg; % blend
  fprintf('[nsd.%s]: random affine transform (tij=[%1.1f,%1.1f], r=%1.2f, sx=%1.2f, sy=%1.2f, kx=%1.1f, ky=%1.1f)\n', mfilename, A_prms);
else
  tmpl = imresize(nsd.util.im2gray(tmplfile), 1);
  img = imresize(nsd.util.im2gray(bgfile), 1);
end


%% Fast hough transform 
[d_img, di_img, fr_img, f_img] = nsd.descriptor(img);
[d_tmpl,di_tmpl,fr_tmpl,f_tmpl] = nsd.descriptor(tmpl);
W = exp(-sqdist(d_img,d_tmpl)/(25));
centroid = [size(tmpl,1)/2,size(tmpl,2)/2]';
dij_tmpl = fr_tmpl(1:2,:)' - repmat(centroid',size(fr_tmpl,2),1);
[ij_max,v_max,Hq] = nsd.hough(W, fr_img(1:2,:)', dij_tmpl, size(img), 8);


%% Display
figure(1); imagesc(Hq, [median(nonzeros(Hq(:))), max(Hq(:))]); title('Hough Voting'); axis image;
figure(2); imagesc(img); colormap(gray); axis image;
hold on; plot(ij_max(1,2),ij_max(1,1),'g+','MarkerSize',30,'LineWidth',5); hold off;


%% Results
if nargin == 0
  fprintf('[nsd.%s]: centroid rmse = %1.3f px\n', mfilename, sqrt(mean((ij_max(1,:)' - (tij+centroid)).^2)));
end
