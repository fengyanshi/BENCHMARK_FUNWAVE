function weights=intWeight(x,type)

% Ensure the vecotr is of the correct length

if isempty(x)
    weights = [];
elseif length(x)==1
    error('must be atleast 2 elements');
else
n = length(x);
delta = x(2)-x(1);
wX = ones(size(x));

if strcmp(type,'simpson')&&mod(n-1,2) ~= 0
    error('length must be odd for simpson rule.');
elseif strcmp(type,'boole')&&mod(n-1,4) ~= 0
    error('length must be = 1 mod 4 for boole rule.');
end

% Generate weighting vector
if strcmp(type,'trap')
    % Trap Rule weighting vector
    wX(2:end-1) = 2;
    weights=delta/2*wX;
elseif strcmp(type,'simpson')
    % Simpsons Rule weighting vector
    wX(2:2:(end-1)) = 4;
    wX(3:2:(end-1)) = 2;
    weights=delta/3*wX;
elseif strcmp(type,'boole')
    % Boole's Rule weighting vector
    wX = 7*wX;
    wX(2:2:(end-1)) = 32;
    wX(3:4:(end-1)) = 12;
    wX(5:4:(end-1)) = 14;
    weights=2*delta/45*wX;
else
    error('type not corrrect');
end
end