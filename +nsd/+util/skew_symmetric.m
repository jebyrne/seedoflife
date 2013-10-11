function [V_hat] = skew_symmetric(V)
%--------------------------------------------------------------------
%
% File: skew_symmetric.m
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
% $Id: skew_symmetric.m 78 2012-07-27 14:11:00Z jebyrne $
%--------------------------------------------------------------------
V_hat = [0 -V(3) V(2); V(3) 0 -V(1); -V(2) V(1) 0];
