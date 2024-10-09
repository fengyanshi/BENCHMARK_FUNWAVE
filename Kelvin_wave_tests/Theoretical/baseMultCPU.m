function z=baseMultCPU(x)
% Perform the matrix vector product for the upper-right/lower-left
% submatrices of the preconditioner matrix
global base;

[n,~] = size(base);
m = length(x);
x = reshape(x,n,m/n);
z = base*x;
z = reshape(z,m,1);