function [ph] = homogenize(p)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne <jebyrne@cis.upenn.edu>
% $Id: homogenize.m 91 2012-12-12 17:03:07Z jebyrne $
%
%--------------------------------------------------------------------------

ph = [p; ones(1,size(p,2))];