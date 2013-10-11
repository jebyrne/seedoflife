function [h] = rectangle(xy, sx, sy, r, linecolor, dotspec)
%--------------------------------------------------------------------------
%
% Copyright (c) 2012 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Inputs
n_pts = size(xy,1);
if ~exist('linecolor','var')
  linecolor = 'green';
end
if ~exist('cornerspec','var')
  cornerspec = 'g.';
end


%% Draw me!
for k=1:n_pts
  A = nsd.util.similarity_transform(xy(k,:),r,1);
  p = [-(sx/2) -(sy/2); -(sx/2) (sy/2); (sx/2) -(sy/2); (sx/2) (sy/2)];
  p = nsd.util.dehomogenize(A*nsd.util.homogenize(p'))';
    
  u = [p(1,:); p(1,:); p(2,:); p(3,:)];
  v = [p(2,:); p(3,:); p(4,:); p(4,:)];
  
  axis ij;  
  hold on; h = line([u(:,1) v(:,1)]', [u(:,2) v(:,2)]','Color',linecolor); hold off;
  hold on; plot(p(:,1),p(:,2),cornerspec); hold off;
  axis equal; 
  drawnow;
end

%   if nargin < 4
%     h = [h line(y(1,:), y(2,:), 'Color', outline_color)];
%   else
%     h = [h fill(y(1,:), y(2,:), fill_color, 'EdgeColor', outline_color)];
%   end
