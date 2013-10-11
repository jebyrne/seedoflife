function [k_inrect] = inrect(ij, rect)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
% $Id: inrect.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

% rect = [xmin ymin width height]
imin = rect(2);
imax = rect(2)+rect(4)-1;
jmin = rect(1);
jmax = rect(1)+rect(3)-1;

% valid ij
%k_valid_i = find((ij(:,1) >= imin) & (ij(:,1) <= imax));
%k_valid_j = find((ij(:,2) >= jmin) & (ij(:,2) <= jmax));
%k_inrect = intersect(k_valid_i, k_valid_j);

% faster
%k_inrect = find((ij(:,1) >= imin) & (ij(:,1) <= imax) & (ij(:,2) >= jmin) & (ij(:,2) <= jmax));
k_inrect = find(all([(ij(:,1) >= imin) (ij(:,1) <= imax) (ij(:,2) >= jmin) (ij(:,2) <= jmax)],2));

