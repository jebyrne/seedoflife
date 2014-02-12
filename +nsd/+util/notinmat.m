function [k_notinmat] = notinmat(matsize, i, j)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
% $Id: notinmat.m 156 2013-04-05 13:59:15Z jebyrne $
%
%--------------------------------------------------------------------------
% Return invalid indices (i,j) for matrix of given size

k_valid = nsd.util.inmat(matsize,i,j);
k_notinmat = setdiff(1:length(i),k_valid);