function [R,T,N,W] = decomposition(H)
%--------------------------------------------------------------------
%
% Authors: Jeffrey Byrne (jbyrne@ssci.com)
%
%--------------------------------------------------------------------


%% 3. Decomposition of the homography matrix
% Use SVD to return sorted singular values 
%[V,D] = eig(H'*H);
[U,S,V] = svd(H'*H);
D = S;

% SO(3) check: footnote 8, pg 138
if (det(V) == -1)
  V = -V;
end

% Use intermediate variables to make this clearer
v1 = V(:,1);
v2 = V(:,2);
v3 = V(:,3);
s1 = D(1,1);
s2 = D(2,2);
s3 = D(3,3);

% MASKS: (5.44) p 136
u1 = ((sqrt(1-s3).*v1) + (sqrt(s1 - 1) .* v3)) ./ (sqrt(s1 - s3));
u2 = ((sqrt(1-s3).*v1) - (sqrt(s1-1).*v3)) ./ sqrt(s1-s3);

% MASKS: p137
U1 = [v2, u1, CrossProdMat(v2)*u1];
U2 = [v2, u2, CrossProdMat(v2)*u2];
W1 = [H*v2, H*u1, CrossProdMat(H*v2)*H*u1];
W2 = [H*v2, H*u2, CrossProdMat(H*v2)*H*u2];

% Solution 1: MASKS table 5.1 (p 137)
solution(1).R = W1*U1'; 
solution(1).N = CrossProdMat(v2)*u1;
solution(1).T = (H-solution(1).R)*solution(1).N;

% Solution 2: MASKS table 5.1 (p 137)
solution(2).R = W2*U2'; 
solution(2).N = CrossProdMat(v2)*u2;
solution(2).T = (H-solution(2).R)*solution(2).N;

% Solution 3: MASKS table 5.1 (p 137)
solution(3).R = solution(1).R;
solution(3).N = - solution(1).N;
solution(3).T = - solution(1).T;

% Solution 4: MASKS table 5.1 (p 137)
solution(4).R = solution(2).R;
solution(4).N = - solution(2).N;
solution(4).T = - solution(2).T;

% Enforce positive depth constraint
valid_solutions = zeros(1,4);
for n=1:4
  % Check solution
  if (solution(n).N'*[0 0 1]' > 0)
    valid_solutions(n) = 1;
  end  
end

% There should be at most two solutions
n = find(valid_solutions == 1);
if (length(n) ~= 2)
  warning(sprintf('[sscv_planar_homography]: Ambiguous solutions: %d', length(n)));
  R = eye(3);
  T = zeros(3,1);
  N = zeros(3,1);
else
  T1 = solution(n(1)).T ./ norm(solution(n(1)).T);
  T2 = solution(n(2)).T ./ norm(solution(n(2)).T);
  %disp(sprintf('[sscv_planar_homography]: Solution 1 (%f,%f,%f), Solution 2 (%f,%f,%f)', T1, T2));
  
  % JB HACK: choose plane most orthogonal to camera 
  if (solution(n(1)).N(3) > solution(n(2)).N(3))
    R = solution(n(1)).R;
    T = solution(n(1)).T;
    N = solution(n(1)).N;    
    T = T ./ norm(T);
  else
    R = solution(n(2)).R;
    T = solution(n(2)).T;
    N = solution(n(2)).N;    
    T = T ./ norm(T);
  end
end

W = DCM2Fixed(R);  % Fixed angles 


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
% $Id: decomposition.m 91 2012-12-12 17:03:07Z jebyrne $
%
%--------------------------------------------------------------------

% Direction cosine matrix --> Fixed angles
b = atan2(-R(3,1),sqrt(R(1,1).^2 + R(2,1).^2));
a = atan2(R(2,1)/cos(b), R(1,1)/cos(b));
g = atan2(R(3,2)/cos(b), R(3,3)/cos(b));

% Output
W = [g b a]';



function [V_hat] = CrossProdMat(V)
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
% $Id: decomposition.m 91 2012-12-12 17:03:07Z jebyrne $
%--------------------------------------------------------------------
V_hat = [0 -V(3) V(2); V(3) 0 -V(1); -V(2) V(1) 0];



