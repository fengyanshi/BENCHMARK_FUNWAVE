% Parameters------------------------------
Fr = 0.7;               % The Froude number
epsilon = 1;            % The pressure strength
freq = 10;              % The sampling frequency
deltaX = 1/freq;        % The mesh spacing
x0 = 0;                 
xEnd = 1000;

x = (x0:deltaX:xEnd);

% The minimum distence to the sensor (must be relatively large to ensure 
% useage of only the single integral term of equation (6.6) is valid,
% y>2/Fr^2 should suffice).
y = 50;

% Calculate solutions
[zeta]=linearZetaWaveTrain(x,y,epsilon,Fr);

save('./Stats/Linear cross-section Fr 0.7.mat','Fr','epsilon','freq','deltaX','x0','x','y','zeta');

%% Computing the spectrogram
% load('./Stats/Linear cross-section Fr 0.7.mat');

% Set up window properties
second = 60;                        % length of teh window in "seconds"
window = 'Harris';                  % Window function
windowLength = round(second*freq);  % Number of data points per window

% Choose the time values for which the spectrogram is computed
ts = (1+ceil(windowLength/2)):100:(length(zeta)-ceil(windowLength/2));

% Compute spectrogram
[tfr]=tfrrsp(zeta,ts,windowLength,window);

%% Organise the data
% Form the time and frequency axes
[fmax,~] = size(tfr);
Nf = floor(fmax/2);
F = freq*(0:(Nf-1))/(2*Nf);
[T,F] = meshgrid(x(ts),F);
P = tfr(1:Nf,:);

% compute the linear dispersion curve
tMax = max(T(:)/y);
t = linspace(sqrt(8),tMax,2001)';
f1 = 1/4*sqrt(2*t.^2+2*t.*sqrt(t.^2-8)+8);
f2 = 1/4*sqrt(2*t.^2-2*t.*sqrt(t.^2-8)+8);

TonY = [t(end:-1:2);t];
omega = [f2(end:-1:2);f1];

%% Plot figure 2(b)
% Ensure matlab properly truncates the colour scheme at the minimum value
cmin = -7;          % Colour axis lower limit
temp = log10(P);
minP = min(min(temp(temp~=-Inf)));
temp(temp<cmin) = cmin;

% Plot specctrogram
figure;
surf(T/y,F*2*pi,temp,'edgecolor','none');
shading interp;
axis tight; 
view(0,90);
yl = ylim;
ylim([yl(1) 5]);
xlabel('t/y','FontSize',13); ylabel('\omega      ','rot',0,'FontSize',15);

% Change the colour scheme from the default
colorbar;
ca = caxis;
caxis([cmin 0]);
% fyshi remove the new setting of colorbar
%cmp = colormap;
%cmp1 = zeros(64,3);
%for i=1:3
%    cmp1(:,i) = interp1(1:64,cmp(:,i),linspace(10,64,64),'spline');
%end
%colormap(cmp1);

% Plot the linear dispersion curve
zMax = get(gca,'zlim');
zMax = zMax(2);
hold on;
plot3(TonY,omega,zMax*ones(size(TonY)),'k');
