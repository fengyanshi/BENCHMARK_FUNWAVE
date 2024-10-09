function mem=calcMem(M,N,band,type)
% Calculate the expended system memory (in Gbs) required to form the
% preconditioner
if strcmp(type,'band')
    kl = (band+1)*(N+1)-1;
    ku = kl;
    n = M*(N+1);
    mem = (2*kl+ku+1)*n+4*M*(N+1)+3*(N+1)^2;
elseif strcmp(type,'dense')
    mem = (M*(N+1))^2+4*M*(N+1);
else 
    mem = (2*M*(N+1))^2;
end
mem = 8*mem/2^30;