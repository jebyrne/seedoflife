function [spyr] = cspyr(img,opt)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne <jebyrne@cis.upenn.edu>
%
%--------------------------------------------------------------------------

%% Options
if ~exist('opt','var')
  opt = nsd.opts().features.spyr;
end


%% Complex steerable pyramid 
spyr.n_bands = opt.n_bands;
spyr.n_scales = opt.n_scales;  
if isempty(spyr.n_scales)
  spyr.n_scales = floor(log2(min(size(img)))) - 2;
end
spyr.opt = opt;
if nsd.util.iscolor(img)
  %imgrey = mat2gray(rgb2gray(img));
  imyuv = nsd.util.imrgb2yuv(img);
  [spyr.pyr,spyr.pind,spyr.steermtx,spyr.harmonics] = buildSCFpyr(double(imyuv(:,:,1)), spyr.n_scales, min(spyr.n_bands-1,15));  % luminance
  for k=2:3
    pyr_chrominance = buildSCFpyr(double(imyuv(:,:,k)), spyr.n_scales, spyr.n_bands-1);
    spyr.pyr = max(spyr.pyr, pyr_chrominance);
  end
else  
% [spyr.pyr,spyr.pind,spyr.steermtx,spyr.harmonics] = buildSCFpyr(double(img), spyr.n_scales, min(spyr.n_bands-1,15)); % steerable basis 
end

sspyr = sepspyr.build(double(img), '9iq', spyr.n_scales, opt.boundary); %warning('sepspyr 9iq')


%% Scales and bands
for z=1:spyr.n_scales
%  [lev,lind] = spyrLev(spyr.pyr,spyr.pind,z);
%  lev = reshape(lev,prod(lind(1,:)),size(spyr.steermtx,2));
  for k=1:spyr.n_bands
    if opt.do_signed_orientation
      r = 2*pi*(k-1)/spyr.n_bands;
    else
      r = pi*(k-1)/spyr.n_bands;
    end    
    %spyr.b{z}{k} = reshape(steer(lev, r, spyr.harmonics, spyr.steermtx), lind(1,:));
    b = sepspyr.steer(sspyr,-r,z); sspyr.br{z}{k} = b{z};  spyr.b{z}{k} = b{z};  % sepspyr, WARNING negative rotation makes big difference?
  end
end
%spyr.lp = pyrLow(spyr.pyr,spyr.pind);
spyr.is_signed = opt.do_signed_orientation;


%% Orientation energy
spyr.n_orientations = spyr.n_bands;
for i=1:spyr.n_scales
  b = [];
  for j=1:spyr.n_bands
    b(:,:,j) = spyr.b{i}{j};
  end
  %b_mag = (2.^(-(i-1))) * abs(b);  warning('testing'); 
  %b_mag = (3.0780.^(-(i-1))) * abs(b);  % common magnitude scale (3.0780)  
  %b_mag = (4.^(-(i-1))) * abs(b);  % common magnitude scale (4)
  b_mag = abs(b); % sepspyr!

%   for j=1:spyr.n_bands
%     figure(1); subplot(1,2,1); imagesc(abs(b_mag(:,:,j))); axis image; colorbar;
%     subplot(1,2,2); imagesc(abs(sspyr.br{i}{j})); axis image; colorbar; %pause();
%   end
    
  b_mag = min(b_mag,opt.oe.max);  % truncate
  b_mag(b_mag < opt.oe.min) = 0;  % minimum activation
  b_ang = angle(b); 
    
  for j=1:spyr.n_bands
    if opt.do_signed_orientation
      %if j > (spyr.n_bands/2), psi0 = 0; else, psi0 = pi; end
      %gamma = cos(b_ang(:,:,j) - psi0).^2;  % phase selection factor [Freeman91]
      %gamma(find(abs(b_ang(:,:,j) - psi0) >= pi/2)) = 0;            
      %gamma = ((b_ang(:,:,j)) ./ (pi));  % linear [-pi,pi] -> [-1,1]
      gamma = 0.5 + (b_ang(:,:,j) ./ (2*pi));  % linear [-pi,pi] -> [0,1]
      %gamma = (1./(1+exp(-b_ang(:,:,j)./(2*pi))));  % sigmoid
      spyr.oe.b{i}{j} = b_mag(:,:,j) .* gamma;
    else
      spyr.oe.b{i}{j} = b_mag(:,:,j);  % orientation energy (cell for fast access)
    end
    %spyr.oe.b{i}{j} = floor(spyr.oe.b{i}{j} .* ((2^opt.n_bits)./opt.oe.max));  % finite precision  
  end

  % Orientation energy visualization
  spyr.oe.img{i} = nsd.util.mat2rgb(b_mag, jet(spyr.n_bands));
end
