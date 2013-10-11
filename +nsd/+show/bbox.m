function [imgout] = bbox(img,bb,opt)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Options
if ~exist('opt','var')
  opt.linewidth = 5; 
  opt.linecolor = [0 1 0];  
end

h = figure(sum(char(mfilename)));
set(h,'NumberTitle','off','Toolbar','none','Menubar','none','Name','Detection');
imagesc(img); axis image; drawnow; 
hold on; rectangle('Position', bb, 'LineWidth',2,'EdgeColor','g'); hold off;
imgout = [];
return;

%% Bounding box image
% % bb = [x y width height]
% n_bb = size(bb,1);
% for k=1:n_bb
%   img_bb = zeros(size(img));
%   img_bb(bb(k,2),bb(k,1):bb(k,3)) = 1;  % top
%   img_bb(bb(k,2):bb(k,4),bb(k,1)) = 1;  % left
%   img_bb(bb(k,4),bb(k,1):bb(k,3)) = 1;  % bottom
%   img_bb(bb(k,2):bb(k,4),bb(k,3)) = 1;  % right
%   img_bb = imdilate(img_bb, ones(opt.linewidth,opt.linewidth));
%   if k == 1
%     img_bb_rgb(:,:,:,k) = ind2rgb(img_bb+1,[0 0 0; opt.linecolor]);
%     img_bb_rgbmask(:,:,:,k) = repmat(img_bb,[1 1 3]);
%   else
%     img_bb_rgb(:,:,:,k) = ind2rgb(img_bb+1,[0 0 0; 1 0 0]) .* (1-img_bb_rgbmask(:,:,:,1));
%     img_bb_rgbmask(:,:,:,k) = repmat(img_bb,[1 1 3]);
%   end
% end
% img_bb_rgb = squeeze(max(img_bb_rgb,[],4));
% img_bb_rgbmask = squeeze(max(img_bb_rgbmask,[],4));
% img_bb_rgbmask = img_bb_rgbmask(1:size(img,1),1:size(img,2),:); % crop
% img_bb_rgb = img_bb_rgb(1:size(img,1),1:size(img,2),:); % crop
% 
% 
% 
% %% Mixture
% imgout = repmat(mat2gray(img), [1 1 3]);  % RGB grayscale
% imgout = imgout .* (1-img_bb_rgbmask) + img_bb_rgb;
% imgout = min(imgout,1);
% 
% 
% %% Show me
% if nargout == 0
%   h = figure(sum(char(mfilename)));
%   set(h,'NumberTitle','off','Toolbar','none','Menubar','none','Name','Detection');
%   imagesc(imgout); axis image; drawnow;
% end
% 
% 
