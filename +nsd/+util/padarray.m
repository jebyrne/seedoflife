function [B] = padarray(A,padsize,padval,paddir)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
% $Id: inrect.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------


switch(padval)
  case {'symmetric','replicate','circular'}
    B = padarray(A,padsize,padval,paddir);
  case {'reflect1'}
    B = padarray(A,padsize+1,'symmetric','both');
    i = padsize(1)+1;
    j = padsize(2)+1;
    u = size(A,1)+1+i;
    v = size(A,2)+1+j;
    B([i u],:) = [];    
    B(:,[j v]) = [];    
    
  otherwise
    if isnumeric(padval)
      B = padarray(A,padsize,padval,paddir);
    else
      error('invalid input');
    end
end
      