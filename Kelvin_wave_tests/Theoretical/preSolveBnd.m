function [z,Flag]=preSolveBnd(y,yscale,Fy,fscale,r)
% Apply the preconditioner matrix to a vector r

global PreSolveCount PreSolveTime data nonLin dataFilepath baseMult;

PreSolveCount = PreSolveCount + 1;
PreSolveTimeTemp = tic;

Flag =0;
l = length(r);
r1 = r(1:l/2);
r2 = r((l/2+1):l);
% z = solveBand(data.J,data.ipvt,r,data.kl);

% even better
y1 = solveBand(data.A,data.ipvtA,r1,1,1,1);
x2 = r2-data.Ccoef*baseMult(y1);
y2 = solveBand(data.D,data.ipvtD,x2,data.kl,data.kl,0);
tempz1 = data.Bcoef*baseMult(y2);
z1 = y1-solveBand(data.A,data.ipvtA,tempz1,1,1,1);
z = [z1;y2];

PreSolveTime = PreSolveTime + toc(PreSolveTimeTemp);
