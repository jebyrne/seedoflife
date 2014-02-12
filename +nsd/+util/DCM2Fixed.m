function [W] = DCM2Fixed(R)
%--------------------------------------------------------------------
%
% File: DCM2Fixed.m
% Authors: Jeffrey Byrne (jbyrne@ssci.com)
%
% Description:  Convert direction cosine matrix to fixed angles. 
% Refer to Craig, "Introduction to Robotics", p47 (2.66)
%
% Inputs:
%   R: direction cosine matrix
%
% Outputs: 
%   W(1)=g: Gamma (roll)
%   W(2)=b: Beta (pitch)
%   W(3)=a: Alpha (yaw)
% 
%   All angle outputs are in radians
%
% $Id: DCM2Fixed.m 78 2012-07-27 14:11:00Z jebyrne $
%
%--------------------------------------------------------------------

% Direction cosine matrix --> Fixed angles
b = atan2(-R(3,1),sqrt(R(1,1).^2 + R(2,1).^2));
a = atan2(R(2,1)/cos(b), R(1,1)/cos(b));
g = atan2(R(3,2)/cos(b), R(3,3)/cos(b));

% Output
W = [g b a]';

