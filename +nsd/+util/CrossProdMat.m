function [V_hat] = CrossProdMat(V)
%--------------------------------------------------------------------
%
% File: CrossProdMat.m
% Authors: Jeffrey Byrne (jbyrne@ssci.com)
%
% Description:  Convert a 3x1 vector V to a 3x3 matrix V_hat such
% that for any 3x1 vector Y: (V x Y) = (V_hat)(Y)
%
% Inputs: 
%   V: 3x1 vector
% 
% Output:
%   V_hat: 3x3 cross product matrix
%
% $Id: CrossProdMat.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------
V_hat = nsd.util.skew_symmetric(V);
