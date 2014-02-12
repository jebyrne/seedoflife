function [spyr] = spyr(img,opt)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne <jebyrne@cis.upenn.edu>
%
%--------------------------------------------------------------------------


%% Steerable pyramid
%filts = 'sp3Filters';  n_orientations = 4;
filts = 'sp5Filters';  n_orientations = 6;
[lo0filt,hi0filt,lofilt,bfilts,steermtx,harmonics] = eval(filts);
z_max = maxPyrHt(size(img),size(lofilt,1));
[spyr.pyr,spyr.pind] = buildSpyr(double(img), z_max, filts);
spyr.n_bands = n_orientations;
spyr.n_scales = z_max;
spyr.filter.lo0filt = lo0filt;
spyr.filter.hi0filt = hi0filt;
spyr.filter.lofilt = lofilt;
spyr.filter.bfilts = bfilts;
spyr.filter.steermtx = steermtx;
spyr.filter.harmonics = harmonics;


%% Bands and scales
for z=1:spyr.n_scales
  for k=1:spyr.n_bands
    spyr.b{z}{k} = spyrBand(spyr.pyr,spyr.pind,z,k);    
  end
end
spyr.lp = pyrLow(spyr.pyr,spyr.pind);


%% Orientation energy 
for i=1:spyr.n_scales
  b = [];
  b(:,:,1) = spyr.b{i}{4};  % orientation bands 0 deg
  b(:,:,2) = spyr.b{i}{5};  % orientation bands 30 deg
  b(:,:,3) = spyr.b{i}{6};  % orientation bands 60 deg
  b(:,:,4) = spyr.b{i}{1};  % orientation bands 90 deg
  b(:,:,5) = spyr.b{i}{2};  % orientation bands 120 deg
  b(:,:,6) = spyr.b{i}{3};  % orientation bands 150 deg
  oe = b.^2;  % orientation energy
  oe = min(oe,opt.oe.max);  % truncate
  oe(oe < opt.oe.min) = 0;  % minimum activation
  for j=1:spyr.n_bands
    spyr.oe.b{i}{j} = oe(:,:,j);  % orientation energy (cell for fast access)
  end

  % Orientation energy visualization
  spyr.oe.img{i} = nsd.util.mat2rgb(oe, opt.oe.cmap);

  % Magnitude 
  spyr.oe.mag{i} = squeeze(sum(b.^2,3));  % over orientation bands
end


%% Display
%showSpyr(spyr.pyr,spyr.pind);
%imagesc(imagesc(spyr.b{1}{1});


