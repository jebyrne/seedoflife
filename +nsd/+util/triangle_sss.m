function [A,B,C] = triangle_sss(a,b,c)
%--------------------------------------------------------------------
%
% Author: Jeffrey Byrne (jbyrne@ssci.com)
% $Id: triangle_sss.m 91 2012-12-12 17:03:07Z jebyrne $
%
%--------------------------------------------------------------------
aa = a;  bb = b; cc = c;
if bb>aa
  aa = bb;
  bb = a;
end
if bb > cc
  mu = cc - (aa-bb);
elseif cc > bb
  mu = bb - (aa-cc);
else
  error('invalid triangle')
end
C = 2*atan( sqrt( ((cc+(aa-bb))*mu)/((aa+(bb+cc))*((aa-cc)+bb)) ) );

% http://www.eecs.berkeley.edu/~wkahan/Triangle.pdf

