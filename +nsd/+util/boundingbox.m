function [bb] = boundingbox(ij)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
% $Id: boundingbox.m 91 2012-12-12 17:03:07Z jebyrne $
%
%--------------------------------------------------------------------------

ijmin = min(ij,[],1);
ijmax = max(ij,[],1);
bb = [ijmin(2) ijmin(1) ijmax(2)-ijmin(2) ijmax(1)-ijmin(1)]; % [x y w h]
