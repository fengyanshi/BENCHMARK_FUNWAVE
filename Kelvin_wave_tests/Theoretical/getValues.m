function [zeta,zetaX,zetaY,phi,phiX,phiY] = getValues(variables,M,N,deltaX,deltaY,x0)
% getValues Returns the surface values. Takes
% the following inputs
% variables - vector of unknowns
% M,N - the number of nodes in the y and x directions respectively
% deltaX,deltaY - the distance between nodes in the x and y directions
% x0 - the smallest x value

% Initialise variables
phi = x0*ones(N,M);
zeta = zeros(N,M);

phiX = reshape(variables(1:(M*(N+1))),(N+1),M);
phi(1,:) = phiX(1,:);
phiX = phiX(2:end,:);
zetaX = reshape(variables((M*(N+1)+1):end),(N+1),M);
zeta(1,:) = zetaX(1,:);
zetaX = zetaX(2:end,:);

% Determine the phi/zeta values
for i=1:N-1
    zeta(i+1,:) = zeta(i,:)+deltaX/2 * (zetaX(i+1,:)+zetaX(i,:));
    phi(i+1,:) = phi(i,:)+deltaX/2 * (phiX(i+1,:)+phiX(i,:));
end
% Determine the phiY/zetaY values

phiY = [zeros(N,1),(phi(:,3:M)-phi(:,1:(M-2)))/(2*deltaY),(phi(:,M)-phi(:,M-1))/(deltaY)];
zetaY = [zeros(N,1),(zeta(:,3:M)-zeta(:,1:(M-2)))/(2*deltaY),(zeta(:,M)-zeta(:,M-1))/(deltaY)];
end