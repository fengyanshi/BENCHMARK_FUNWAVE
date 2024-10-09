function [zeta]=linearZetaWaveTrain(x,y,epsilon,Fr)
% linearZetaWaveTrain Computes the surface height of the linear wave train
% for flow past a Gaussian pressure distribution over a rectangular mesh.
% The wave train is found by computing the single integral term of equation
% (6.6) of the thesis
% x - x values of the mesh
% y - y values of the mesh
% epsilon - strength of teh pressure
% Fr - the Froude number

N = length(x);
M = length(y);

% For all y values
for j=1:M
    
    % Determine the required grid spacing for evaluating the integral
    n = 2;
    acc = 2;
    p = acos(nthroot(1/(acc*4*pi^2*Fr^4*log(10)),4));
    temp = pi*(pi-2*p)^2-4*n*max(y(j)/Fr^2,1);
    dp = -(pi-2*p)/2*(temp+2*sqrt(-n*max(y(j)/Fr^2,1)*temp))/temp;
    dpx = -(1/2)*(-2*p+pi)^2*pi/(pi^2-2*pi*p-2*n*max(x(end)/Fr^2,1));
    psiDet = round(max(ceil(pi/dp),ceil(pi/dpx)));
    
    % Generate the grid and the weightings for integration
    psi = linspace(-pi/2,pi/2,psiDet)';
    psiWeight = intWeight(psi,'trap');

    % Store often used values
    cosPsi = cos(psi);
    cosPsi([1,psiDet]) = 0;
    sinPsi = sin(psi);
    k0 = 1./(cosPsi.^2);
    k0Sq = k0.^2;
    
    % The fourier transform of the pressure distribution

% fyshi: this P is p(k,theta), the difference is F_L, not Fr. in fact there s a typo in (11), 
% x^2 should be (x^2 + y^2)
% -------------------

    P = 1/(Fr^4*pi)*exp(-k0Sq/(4*Fr^4*pi^2));
    parfor k=1:N
        % If statment represent the Heaviside function
        if x(k)<0
            zeta(k,j) = 0;
        else
            % Compute integral
            temp = k0Sq.*P.*sin((cosPsi*x(k)-sinPsi*y(j)).*k0);
            temp(isnan(temp)) = 0;
            zeta(k,j) = -epsilon/pi*real(sum(psiWeight.*temp));
        end
    end
end



end