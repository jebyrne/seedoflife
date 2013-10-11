function [p] = dehomogenize(ph)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne <jebyrne@cis.upenn.edu>
% $Id: dehomogenize.m 91 2012-12-12 17:03:07Z jebyrne $
%
%--------------------------------------------------------------------------

d = size(ph,1)-1;
p = ph(1:end-1,:) ./ repmat(ph(end,:),d,1);
