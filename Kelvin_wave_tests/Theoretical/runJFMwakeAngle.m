% Parameters------------------------------
Fr = 0.9;               % The Froude number
singType = false;       % The type of singularity, true = source, false=dipole
N = 735;                % The mesh size
M = 237;                %

% list of epsilons for which the wake angle is measured
AllEps = [0.02:0.06:0.92];

% The folder that the surface data is saved to  
folder = './Stats/Dipole-Fr09/';
% Calculate solutions
[peaks,gofs] = measureAllPeaks(folder,N,M,Fr,AllEps,singType);

figure;plot(AllEps,peaks,'o');
