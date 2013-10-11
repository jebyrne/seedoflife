function [k_inbbox] = inbbox(ij, bbox)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
% $Id: inrect.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

% bbox = [xmin ymin xmax ymax] = [jmin imin jmax imax]
imin = bbox(2);
imax = bbox(4);
jmin = bbox(1);
jmax = bbox(3);

k_inbbox = find(all([(ij(:,1) >= imin) (ij(:,1) <= imax) (ij(:,2) >= jmin) (ij(:,2) <= jmax)],2));

