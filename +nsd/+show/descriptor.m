function [] = descriptor(img,di,h,linewidth)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Inputs
if ~exist('linewidth','var')
  linewidth = 1;
end


%% Figure setup
if ~exist('h','var')
  h = figure(sum(char(mfilename)));
  set(h,'NumberTitle','off','Toolbar','none','Menubar','none','Name','Descriptors');
end


%% Show me
imagesc(img); colormap(gray); axis image; axis off;
hold on;
n_desc = size(di.fr,2);
for i=1:n_desc
  % Nested support
  for j=1:di.n_scales
    xy = [di.fr(2,i) di.fr(1,i)];
    r = di.fr(4,i);
    s = di.r_nesting(j);
    h = nsd.show.circle(xy',s,'b',[],256);
    set(h,'LineWidth',linewidth);
    
    % Shape support
%     for k=1:di.n_shapesupport
%       mu = nsd.util.ij2xy(di.dij_shapesupport{j}(k,:) + di.fr(1:2,i)',2)';
%       if j > 2
%         r = ((di.m_scaledsupport(j-1))/2);
%         h=draw_circle(mu,r,'b');
%         set(h,'LineWidth',linewidth);
%       end
%     end
  end  
end
hold off;
drawnow;



%% ARCHIVE
% Quiver
% [u,v] = pol2cart(fr(4,:)',r_scale*fr(3,:)');
% hold on; quiver(fr(2,:),fr(1,:),u(1:1:end)',v(1:1:end)',0,'y-'); hold off;
% axis equal;
% drawnow;

%fr_vlfeat = [fr(2,:); fr(1,:); fr(3,:); fr(4,:)];  % ij -> xy for vlfeat
%hold on;  h1 = vl_plotframe(fr_vlfeat); hold off;
%set(h1,'color','y','linewidth',3);

