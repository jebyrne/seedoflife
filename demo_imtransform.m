function [] = demo_imtransform(imfile)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne
%
%--------------------------------------------------------------------------
close all;


%% Inputs
if nargin == 0
  imfile = 'peppers.png';
end
im = imread(imfile); 
if nsd.util.iscolor(im)
  im = mat2gray(rgb2gray(im));
end


%% Random transformation
txy = 16*(randn(2,1));  % translation in (x,y)
r = sign(randn)*rand*(pi/4);  % rotation [0,2pi]
sx = randn*0.5 + 1-0.125;  % [0.5,1]
sy = randn*0.5 + 1-0.125;  % [0.5,1]
kx = 0.125*randn + nsd.util.randsign()*0; % [-0.5,0.5]
ky = 0.125*randn + nsd.util.randsign()*0; % [-0.5,0.5]


%% Affine transformation
[A,xx,prms] = nsd.util.affine_transform(txy,r,sy,sx,kx,ky);  
[im_xform,B] = nsd.util.imtransform(im,A);
xyh = [rand(1,16)*size(im,2); rand(1,16)*size(im,1); ones(1,16)];  % homogeneous point in im
xyh_xform = B*xyh;  % corresponding point in im_xform

figure(1); subplot(3,2,1); imagesc(im); colormap(gray); axis image;
hold on; plot(xyh(1,:),xyh(2,:),'r.','MarkerSize',20); hold off;
figure(1); subplot(3,2,2); imagesc(im_xform); colormap(gray); axis image;
hold on; plot(xyh_xform(1,:),xyh_xform(2,:),'r.','MarkerSize',20); hold off;
title('affine transform');
fprintf('[nsd.%s]: random affine transform - tij=[%1.1f,%1.1f], r=%1.2f, s=[%1.2f,%1.2f], k=[%1.2f,%1.2f]\n', mfilename, prms.txy(1),prms.txy(2),prms.r, prms.sx,prms.sy,prms.kx,prms.ky);


%% Similarity transform
[A,xx,prms] = nsd.util.similarity_transform(txy,r,min(sx,sy)); 
[im_xform,B] = nsd.util.imtransform(im,A);
xyh = [rand(1,16)*size(im,2); rand(1,16)*size(im,1); ones(1,16)];  % homogeneous point in im
xyh_xform = B*xyh;  % corresponding point in im_xform

figure(1); subplot(3,2,3); imagesc(im); colormap(gray); axis image;
hold on; plot(xyh(1,:),xyh(2,:),'g.','MarkerSize',20); hold off;
figure(1); subplot(3,2,4); imagesc(im_xform); colormap(gray); axis image;
hold on; plot(xyh_xform(1,:),xyh_xform(2,:),'g.','MarkerSize',20); hold off;
title('similarity transform');
fprintf('[nsd.%s]: random similarity transform - tij=[%1.1f,%1.1f], r=%1.2f, s=[%1.2f]\n', mfilename, prms.txy(1),prms.txy(2),prms.r, prms.s);


%% Homography
[H] = nsd.homography.random(size(im,1),size(im,2),4);
[im_xform,B] = nsd.util.imtransform(im,H);
xyh = [rand(1,16)*size(im,2); rand(1,16)*size(im,1); ones(1,16)];  % homogeneous point in im
xyh_xform = B*xyh;  % corresponding point in im_xform
xy_xform = [xyh_xform(1,:)./xyh_xform(3,:); xyh_xform(2,:)./xyh_xform(3,:)];  % dehomogenize

figure(1); subplot(3,2,5); imagesc(im); colormap(gray); axis image;
hold on; plot(xyh(1,:),xyh(2,:),'b.','MarkerSize',20); hold off;
figure(1); subplot(3,2,6); imagesc(im_xform); colormap(gray); axis image;
hold on; plot(xy_xform(1,:),xy_xform(2,:),'b.','MarkerSize',20); hold off;
title('projective transform');
fprintf('[nsd.%s]: random homography\n', mfilename);
H
