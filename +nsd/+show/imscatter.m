function [imgout] = imscatter(img,ij,w,h,s)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Figure setup
if ~exist('h','var')
  h=figure(sum(char(mfilename)));
  set(h,'NumberTitle','off','Toolbar','none','Menubar','none','Name','Scatterplot');
end
if ~exist('w','var') || isempty(w)
  w = ones(size(ij,1),1);
else
  w = mat2gray(w); % [0,1]
end
if ~exist('s','var')
  s = 10;
end



%% Visualization
if ~isempty(img)
  [k,k_valid] = nsd.util.sub2ind(size(img),ij(:,1),ij(:,2));  % do not include invalid ij
else
  k_valid = 1:size(ij,1);
end
%[xx, k_sort] = sort(k);  
%ij = ij(k_valid(k_sort),:);
ij = ij(k_valid,:);  % do not sort
w = w(k_valid); 

n_points = size(ij,1);
S = round(s*w)+1;
C = jet(n_points);


%% Image
if ~isempty(img)
  figure(h); imagesc(img); axis image; colormap(gray); axis off; 
end
hold on; 
scatter(ij(:,2),ij(:,1),S,C,'filled'); axis ij; 
hold off;
drawnow;

