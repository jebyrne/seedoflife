function [imgout] = detection(img,ij,ij_true)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------


%% Green "+" image 
markerwidth = 2;
markerheight = 12;
imgout = repmat(img,[1 1 3]);  % grayscale
if exist('ij_true','var')
  ij_true = round(ij_true);
end
try
  ij = round(ij);
  imgout(ij(1)-markerwidth:ij(1)+markerwidth,ij(2)-markerheight:ij(2)+markerheight,1) = 0;
  imgout(ij(1)-markerwidth:ij(1)+markerwidth,ij(2)-markerheight:ij(2)+markerheight,2) = 1;
  imgout(ij(1)-markerwidth:ij(1)+markerwidth,ij(2)-markerheight:ij(2)+markerheight,3) = 0;
  imgout(ij(1)-markerheight:ij(1)+markerheight,ij(2)-markerwidth:ij(2)+markerwidth,1) = 0;
  imgout(ij(1)-markerheight:ij(1)+markerheight,ij(2)-markerwidth:ij(2)+markerwidth,2) = 1;
  imgout(ij(1)-markerheight:ij(1)+markerheight,ij(2)-markerwidth:ij(2)+markerwidth,3) = 0;
catch ME
  warning('centroid too close to border for this hacky display');
end


%% Show me
if nargout == 0
  h = figure(sum(char(mfilename)));
  set(h,'NumberTitle','off','Toolbar','none','Menubar','none','Name','Detection');
  imagesc(img); colormap(gray); axis image; axis off;
  if exist('ij_true','var')
    hold on; plot(ij_true(2), ij_true(1), 'r+', 'MarkerSize', 20, 'LineWidth',5); hold off;
  end
  hold on; plot(ij(2), ij(1), 'g+', 'MarkerSize', 20, 'LineWidth',5); hold off;
  drawnow;
end



