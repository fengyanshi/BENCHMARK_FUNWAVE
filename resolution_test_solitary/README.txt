
folders
1)
/dx_p1_ht_p01/ dx=0.1, ht = 0.01*ht0, Cd = 0.0
2)
/dx_p1_ht_p1/ dx=0.1, ht = 0.1*ht0, Cd = 0.0
3)
/dx_p1_ht_p3/ dx=0.1, ht = 0.3*ht0, Cd = 0.0
4)
/dx_p05_ht_p3/ dx=0.05, ht = 0.3*ht0, Cd = 0.0
5)
/dx_p1_ht_p1/ dx=0.05, ht = 0.1*ht0, Cd = 0.0
6) 
/dx_p1_ht_p01_frc/ dx=0.1, ht = 0.01*ht0, Cd = 0.01
7)
/dx_p1_ht_p1_frc/ dx=0.1, ht = 0.1*ht0, Cd = 0.01

Use /postprocessing/*.m to draw figures

Conclusions
1) As ratio = Ht/Depth becomes large, tails are generated due to inconsistency between FUNWAVE equations and the weakly nonlinear solution 

2) the grid should resolve the solitary peak, at least 20 points over 80% peak region, or 80 points over the whole length. So case 3) is not proper model setup

3) When bottom friciton applied, wave is damped significantly. See cases 6) and 7)