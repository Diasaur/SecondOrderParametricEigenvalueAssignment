function [A,B,E] = SecOrd2Descriptor(A_2,A_1,A_0,B_0)
% Takes a second order system (of the form A2 dd + A1 dx + A0 x = u) and 
% outputs an LTI system d[x;dx]=A[x;dx]+Bu

temp = size(A_2);
n = temp(1);  % original dimension
temp = size(B_0);
m = temp(2);

% calculating system matrices (A,B)
A = [zeros(n), eye(n); -A_0, -A_1];  % system dynamics matrix
B = [zeros(n,m); B_0];  % input matrix
E = [eye(n), zeros(n); zeros(n), A_2];

end