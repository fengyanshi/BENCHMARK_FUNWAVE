function setTo0(x)
global PreSolveCount PreSolveTime SawPreRunCount SawPreRunTime...
    FDJacRunCount FDJacRunTime LinAlgJacRunCount LinAlgJacRunTime...
    LinAlgPreRunCount LinAlgPreRunTime LinFDPreRunCount LinFDPreRunTime...
    LinSawPreRunCount LinSawPreRunTime NonLinFuncRunCount NonLinFuncRunTime...
    SawFuncRunCount SawFuncRunTime LinFuncRunCount LinFuncRunTime LinPreRunTime...
    oldPreSolveCount PreSetupCount;

PreSetupCount = 0;
PreSolveCount = 0;
PreSolveTime = 0;
SawPreRunCount = 0;
SawPreRunTime = 0;
FDJacRunCount = 0;
FDJacRunTime = 0;
LinAlgJacRunCount = 0;
LinAlgJacRunTime = 0;
LinAlgPreRunCount = 0;
LinAlgPreRunTime = 0;
LinFDPreRunCount = 0;
LinFDPreRunTime = 0;
LinSawPreRunCount = 0;
LinSawPreRunTime = 0;
NonLinFuncRunCount = 0;
NonLinFuncRunTime = 0;
SawFuncRunCount = 0;
SawFuncRunTime = 0;
LinFuncRunCount = 0;
LinFuncRunTime = 0;
if x    
    LinPreRunTime = 0;
end
oldPreSolveCount = 0;
end