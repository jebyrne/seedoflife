function [h] = matching(imobs,imref,ij_ref,ij_refinobs,ij_obs,h,w,smax)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Figure setup
if ~exist('h','var') || isempty(h)
  %h=figure(sum(char(mfilename)));
  h=figure;
  set(h,'NumberTitle','on','Name','Matching');
end
if ~exist('w','var') || isempty(w)
  w = ones(size(ij_ref,1),1);
else
  w = mat2gray(w); % [0,1]
end
if ~exist('smax','var') || isempty(smax)
  smax = 10;
end
if ~exist('do_pause','var') || isempty(do_pause)
  do_pause = 0;
end


%% Visualization
[k_img,k_valid] = nsd.util.sub2ind(size(imobs),ij_refinobs(:,1),ij_refinobs(:,2));  % do not include invalid ij not in observation
%k_sort = 1:length(k_img); % colors in reference order
[xx, k_sort] = sort(k_img);  % colors in observed order
ij_refinobs = ij_refinobs(k_valid(k_sort),:);
ij_ref = ij_ref(k_valid(k_sort),:); 
w = w(k_valid);

n_points = size(ij_ref,1);
S = round(smax*w)+5;
C = jet(n_points);


%% Observation
figure(h);
subplot(1,2,1); 
imagesc(mat2gray(imobs)); axis image; colormap(gray); axis off; title('Observation');
hold on; 
if ~isempty(ij_obs)
  scatter(ij_obs(:,2),ij_obs(:,1),2,'k','filled'); axis ij; 
end
scatter(ij_refinobs(:,2),ij_refinobs(:,1),S,C,'filled'); axis ij; 
hold off;


%% Reference
subplot(1,2,2);  
imagesc(mat2gray(imref)); axis image; colormap(gray); axis off; title('Reference');
hold on;
scatter(ij_ref(:,2),ij_ref(:,1),smax*ones(size(S)),C,'filled'); axis ij;   
hold off;


%% Done
drawnow;

