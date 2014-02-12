function [imgout] = assignment(img,tmpl,ij_img,ij_tmpl)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------
imgout = nsd.show.matching(img,tmpl,ij_img,ij_tmpl);

