function [peaks,gofs] = measureAllPeaks(folder,N,M,Fr,AllEps,singType)
% measureAllPeaks The function that returns the apparent wake angle (peaks)
% and goodness of fit (gof) for the line of best fit for a set of different
% source strngths epsilon
% folder - The folder were the saved surface data is stored
% M,N - the number of nodes in the y and x directions respectively
% Fr - the Froude number
% AllEps - a vector of epsilons
% singType - (bool) The type of singularity, true=source, false=dipole

n = length(AllEps);
if singType
    singularity = 'Source';
    singPara = 'eps';
else
    singularity = 'Dipole';
    singPara = 'mu';
end
peaks = zeros(n,1);
gofs = cell(n,1);
for k=1:n
    
    % look in the given folder for the correct surface files
    listings = dir([folder,singularity,' ',num2str(N),'x',num2str(M),' ',singPara,' ',num2str(AllEps(k)),' F ',num2str(Fr),' *']);
        
    if ~isempty(listings)
        % Check if more than one surface file has the same nondimensional
        % parameters, if so load the latest solution
        ll = length(listings);
        if ll==1
            load([folder,listings(1).name]);
        elseif ll>1
            dates = zeros(ll,1);
            for j=1:ll
                dates(j)=listings(j).datenum;
            end
            [~,ind] = max(dates);
            load([folder,listings(ind).name]);
        end
        
        % Find the apparent wake angle and teh goodness of fit for the line
        % of best fit
        [pA,gof]=getPeakAngle(yDash,M,N,deltaX,deltaY,x0,Fr);
        peaks(k) = pA;
        gofs{k} = gof;
    else
        peaks(k) = NaN;
    end
end