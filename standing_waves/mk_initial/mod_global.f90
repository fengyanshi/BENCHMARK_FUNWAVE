     Module GLOBAL

     implicit none

! fyshi added save 04/01/2013
     SAVE

     ! define precision

     integer, parameter :: SP=8


     ! define parameters
     real(SP), parameter :: pi=3.141592653
     real(SP), parameter :: Small=0.000001
     real(SP), parameter :: Large=10000000.0
     real(SP), parameter :: Grav=9.81
     real(SP), parameter :: Zero=0.0
     real(SP), parameter :: One=1.0
     real(SP), parameter :: Rho0=1000.0
     real(SP), parameter :: RhoA=1.20
     real(SP), parameter :: Kappa=0.41

! fyshi change to integer 12/15/2011
     integer,  parameter :: MaxNumFreq=100
     integer,  parameter :: MaxNumDir=100

     ! ghost cells (>=1)
     integer, parameter :: Nghost=2

     ! define characters
     character(len=80) :: TITLE
     character(len=80) :: RESULT_FOLDER
     character(len=80) :: HIGH_ORDER
     character(len=80) :: TIME_ORDER
     character(len=80) :: WaveMaker
! fyshi add character for boundary condition 11/02/2012
     character(len=80) :: BOUNDARY
     character(len=80) :: DEPTH_TYPE
     character(len=80) :: dt_constraint
     character(len=80) :: CONVECTION








! fyshi create temporary array
     real(SP), dimension(:,:), allocatable :: tmp_2d_1,tmp_2d_2







     ! define output logical parameters
     logical :: ANA_BATHY,NON_HYDRO,VISCOUS_FLOW,SPONGE_ON,OUT_H,OUT_E,OUT_U,OUT_V,OUT_W,OUT_P, &
                OUT_K,OUT_D,OUT_S,OUT_C,OUT_B,OUT_A,OUT_F,OUT_T,OUT_G,OUT_I,PERIODIC_X,PERIODIC_Y, &
                WAVE_AVERAGE_ON,ADV_HLLC,BAROTROPIC,RIGID_LID,BED_CHANGE,EXTERNAL_FORCING,STATIONARY

     ! variables
     integer :: It_Order,Ibeg,Iend,Iend1,Jbeg,Jend,Jend1,Kbeg,Kend,Kend1,PX,PY,IVturb,IHturb,  &
                Mglob,Nglob,Kglob,Mloc,Nloc,Kloc,Mloc1,Nloc1,Kloc1,Icount,RUN_STEP,Ivgrd,SIM_STEPS,Ibot, &
                NumFreq,NumDir,NSTAT,WaveheightID


     integer :: Bc_X0,Bc_Xn,Bc_Y0,Bc_Yn,Bc_Z0,Bc_Zn
     real(SP) :: dt,dt_old,dt_min,dt_max,dt_ini,dx,dy,Theta,CFL,  &
                VISCOUS_NUMBER,MinDep,TIME,TOTAL_TIME,Plot_Intv,  &
                 Screen_Intv,Screen_Count,Plot_Count,Visc,Cvs,Chs,Zob,Tke_min,  &
                 Eps_min,Cmut_min,Cd0,Plot_Start,Plot_Intv_Stat, &
                 Plot_Count_Stat,xstat(20),ystat(20),Wave_Ave_Start,Wave_Ave_End,Schmidt,TRamp,Grd_R
     real(SP) :: Amp_Wave,Per_Wave,Dep_Wave,Theta_Wave,Freq(MaxNumFreq),  &
                Dire(MaxNumDir),Wave_Spc2d(MaxNumDir,MaxNumFreq), &
                 Random_Phs(MaxNumDir,MaxNumFreq),Hm0,Tp,Freq_Min,  &
                Freq_Max,Jon_Spc(MaxNumFreq),RanPhs(MaxNumFreq)
     real(SP) :: Sponge_West_Width,Sponge_East_Width,Sponge_South_Width,  &
                Sponge_North_Width,R_Sponge,A_Sponge, &
                 Xsource_West,Xsource_East,Ysource_Suth,Ysource_Nrth
     real(SP), dimension(3) :: ALPHA,BETA

! fyshi added time series boundary condition 12/17/2011
       CHARACTER(LEN=80) :: BoundaryFile,WHAT
       INTEGER :: NumTimeData
       INTEGER :: icount_tide = 1
       REAL(SP),DIMENSION(:),ALLOCATABLE :: DataU_L,DataEta_L,DataSal_L,DataTem_L
       REAL(SP),DIMENSION(:),ALLOCATABLE :: DataU_R,DataEta_R,DataSal_R,DataTem_R
       REAL(SP),DIMENSION(:),ALLOCATABLE :: TimeData
       REAL(SP),DIMENSION(:),ALLOCATABLE :: Z_pct_West,Z_pct_East

! fyshi added bathymetry file 04/13/2012
       CHARACTER(LEN=80) :: Depth_File

     ! real arrays
     real(SP), dimension(:), allocatable :: x,xc,y,yc,sig,  &
                dsig,sigc,Ein_X0,Din_X0,Ein_Xn,Din_Xn, &
                Ein_Y0,Din_Y0,Ein_Yn,Din_Yn
     real(SP), dimension(:,:), allocatable :: Ho,H,Hc,HCG,Hc0,Hfx,Hfy,  &
                DeltH,DeltHo,Delt2H,DelxH,DelyH,D,D0,  &
                Eta,Eta0,Eta00, &
                SourceX,SourceY,SourceC,DxL,  &
                DxR,DyL,DyR,EtaxL,EtaxR,EtayL,EtayR, &
                DelxEta,DelyEta,DelxD,DelyD,Uin_X0,  &
                Vin_X0,Win_X0,Uin_Xn,Vin_Xn, &
                Win_Xn,Bc_Prs,Sponge,Setup,WaveHeight,  &
                Uin_Y0,  &
                Vin_Y0,Win_Y0,Uin_Yn,Vin_Yn, &
                Win_Yn, &
                Umean,Vmean,Emax,Emin
     real(SP), dimension(:,:,:), allocatable :: U,V,W,U0,V0,W0,  &
                U00,V00,W00,Omega,P,DU,DV,DW,DU0,DV0,DW0, &
                UxL,UxR,VxL,VxR,WxL,WxR,DUxL,DUxR,DVxL,DVxR,DWxL, &
                DWxR,UyL,UyR,VyL,VyR,WyL,WyR,DUyL,DUyR,DVyL,DVyR,DWyL,DWyR, &
                UzL,UzR,VzL,VzR,WzL,WzR,OzL,OzR,SxL,SxR,SxS,SyL,SyR,SyS,ExL,ExR,FxL, &
                FxR,GxL,GxR,HxL,HxR,EyL,EyR,FyL,FyR,GyL,GyR,HyL,HyR,Ex,Ey,Fx, &
                Fy,Fz,Gx,Gy,Gz,Hx,Hy,Hz,DelxU,DelyU,DelzU,DelxV,DelyV,DelzV, &
                DelxW,DelyW,DelzW,DelxDU,DelyDU,DelxDV,DelyDV,DelxDW,DelyDW, &
                DelzO,Uf,Vf,Wf,Cmu,Cmuht,Cmuvt,Diffxx,Diffxy,Diffxz,Diffyx,  &
                Diffyy,Diffyz,Diffzx,Diffzy,Diffzz,DelxSc,DelySc,Rho,Rmean,Tke,Eps,Skl, &
                DTke,DEps,DTke0,DEps0,Prod_s,Prod_b,Lag_Umean,Lag_Vmean,Lag_Wmean, &
                Euler_Umean,Euler_Vmean,Euler_Wmean,DRhoX,DRhoY,ExtForceX,ExtForceY, &
                                                UpWp

     ! integer arrays
     integer, dimension(:,:), allocatable :: Mask,Mask_Struct,Mask9,Brks,Num_Zero_Up
     
     ! poisson solvers
     integer  :: itmax,isolver,neqns
     real(SP) :: tol
     real(SP), dimension(:),   allocatable :: Rhs
     integer,  dimension(:),   allocatable :: JCoef
     real(SP), dimension(:,:), allocatable :: Coef

! fyshi add initial sali and temp conditions
     CHARACTER(LEN=80) :: INI_SALI_INPUT,INI_TEMP_INPUT
     CHARACTER(LEN=80) :: INI_SALI_FILE,INI_TEMP_FILE     
     REAL(SP) :: INI_SALI=35.0
     REAL(SP) :: INI_TEMP=0.0

! fyshi add tidal current low pass
     LOGICAL :: TID_LOW_PASS=.FALSE.

! fyshi added nesting option 05/15/2013



     End Module GLOBAL

