!------------------------------------------------------------------------------------
!
!      FILE main.F
!
!      This file is part of the FUNWAVE-TVD program under the Simplified BSD license
!
!-------------------------------------------------------------------------------------
! 
!    Copyright (c) 2016, FUNWAVE Development Team
!
!    (See http://www.udel.edu/kirby/programs/funwave/funwave.html
!     for Development Team membership)
!
!    All rights reserved.
!
!    FUNWAVE_TVD is free software: you can redistribute it and/or modify
!    it under the terms of the Simplified BSD License as released by
!    the Berkeley Software Distribution (BSD).
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions are met:
!
!    1. Redistributions of source code must retain the above copyright notice, this
!       list of conditions and the following disclaimer.
!    2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!  
!    The views and conclusions contained in the software and documentation are those
!    of the authors and should not be interpreted as representing official policies,
!    either expressed or implied, of the FreeBSD Project.
!  
!-------------------------------------------------------------------------------------
Program FUNWAVE_TVD
!-------------------------------------------------------------------------------------
!    VERSION 3.0
!
!    MAIN - READ_INPUT
!         - INDEX
!         - ALLOCATE_VARIABLES
!         - INITIALIZATION 
!         DO LOOP
!             - VARIABLE_UPDATING
!             - EXCHANGE
!             - ESTIMATE_DT
!             - RUNGE-KUTTA
!             - DISPERSION
!             - FLUXES          
!             - SourceTerms      
!             - ESTIMATE_HUV
!             - EXCHANGE 
!             - FLUXES                 
!             - SourceTerms again
!
!             - STATISTICS
!             - PREVIEW
!          ENDDO LOOP
!-------------------------------------------------------------------------------------
! ** OPEN FILES **
!  (1): read input, (2): output, (3): log, (4): !write/read hotstart
!-------------------------------------------------------------------------------------
! ** HOT START DATA **
!   NOTE: read input.txt first, if HOT_START, then read  
        ! -- dimension
! Mloc,Nloc,Mloc1,Nloc1
! Nghost
! Ibeg,Iend,Jbeg,Jend,Iend1,Jend1
!   NOTE: need to confirm if the saved data is consistent with input.txt
        ! -- time
! TIME
! TOTAL_TIME
! PLOT_INTV
! PLOT_COUNT
! SCREEN_INTV
! SCREEN_COUNT
! HOTSTART_INTV
! ICOUNT
        ! spacing
! DX,DY
        ! -- physics
! DISPERSION
! Gamma1
! a1,a2,b1,b2
! SWE_ETA_DEP
        ! -- numerics
! Time_Scheme
! HIGH_ORDER
! CONSTR
! CFL
! FroudeCap
! DISP_TIME_LEFT
        ! -- wet-dry
! MinDepth,MinDepthfrc

        ! -- depth
! DEPTH
! DEPTHx
! DEPTHy
        ! variables
! U
! V
! if (.NOT.DISP_TIME_LEFT)THEN
! U0
! V0
! endif
! Ubar
! Vbar
! ETA 
! H
! MASK
! MASK9
! MAST_STRUC
!
       ! -- wavemaker
! if (WAVEMAKER is WK_IRR)
! turns out the data for Cm Sm too large, calculate it when hotstart
!
! if (WAVEMAKER is WK_REG)
! D_gen
! Beta_gen
! rlamda
! 
!
!-------------------------------------------------------------------------------------
     USE GLOBAL
     






     USE TIDE_MODULE
     USE TIME_SPECTRA_MODULE

     IMPLICIT NONE

!     INTEGER::ISTAGE ! moved to mod_global 09/12/2017
!     REAL(SP) :: tbegin,tend  ! moved to mod_global 07/29/2016


     CALL MPI_INIT ( ier )




     CALL READ_INPUT



	![ykchoi(14.12.24.)
	IF(INI_UVZ)THEN
        TIME=HotStartTime
      ENDIF
	!ykchoi(14.12.24.)]

     CALL INDEX



! allocate variables
     CALL ALLOCATE_VARIABLES



     CALL INITIALIZATION

! time dependent wave spectra

     IF(WaveMaker(1:12)=='TIME_SPECTRA') THEN
       CALL TIME_SPECTRA_INITIAL
     ENDIF

     CALL TIDE_INITIAL

     IF(WaveMaker(1:12)=='TIME_SPECTRA') THEN
!    time spectra only use tide abs_gen and tide data type
     CALL TIME_SPECTRA_PROTECTION
     ENDIF
      


















! time integration

     ! record wall time

     if(myid == 0) tbegin = MPI_Wtime( )



   DO WHILE (TIME<TOTAL_TIME)

!     move output here to get the initial condition 11/27/2018

      CALL OUTPUT

     IF(WaveMaker(1:7)=='LEF_SOL')THEN
       CALL SOLITARY_WAVE_LEFT_BOUNDARY
     ENDIF   

! update three variables
     Eta0=Eta
     Ubar0=Ubar
     Vbar0=Vbar  


! previous update_mask
   



     CALL EXCHANGE



 
     CALL ESTIMATE_DT(Mloc,Nloc,DX,DY,U,V,H,MinDepthFrc,DT,CFL,TIME)


  
! U0, V0 are moved to following part due to computation of Ut, Vt.
	U0=U   !ykchoi(15. 08. 06.)
	V0=V   !ykchoi







       IF(TideBcType(1:4)=='DATA')THEN
         CALL TIDE_DATA
       ENDIF



       IF(WaveMaker(1:12) == 'TIME_SPECTRA') THEN
         CALL TIME_SPECTRA_INTERPOLATION
       ENDIF






     ! 3-ORDER RUNGE-KUTTA TIME STEPPING
     DO ISTAGE=1,3

       IF(DISPERSION)THEN
         CALL Cal_Dispersion
       ENDIF 



       CALL FLUXES



       CALL SourceTerms   ! put sourceterms after fluxes in order to get eta_t



       CALL ESTIMATE_HUV(ISTAGE) 

! etascreen was added, update_mask was moved here from outside RK 
       CALL UPDATE_MASK

       IF(TIDAL_BC_ABS)THEN
         CALL TIDE_BC
       ENDIF





       CALL WAVE_BREAKING




       
       CALL EXCHANGE


       

         

       IF(WaveMaker(1:3)=='ABS' .OR. WaveMaker(1:12)=='TIME_SPECTRA') THEN
         CALL ABSORBING_GENERATING_BC



       ENDIF

       IF(WaveMaker(1:11)=='LEFT_BC_IRR') THEN
         CALL IRREGULAR_LEFT_BC
       ENDIF

       IF(DIRECT_SPONGE)THEN
           CALL SPONGE_DAMPING
       ENDIF

     ENDDO



     CALL MIXING_STUFF






!  find maximum eta velocity 

      IF (OUT_Hmax.OR.OUT_Hmin.OR.OUT_Umax.OR.OUT_MFmax.OR.OUT_VORmax.OR.OUT_Time)THEN
        CALL MAX_MIN_PROPERTY
      ENDIF        

      CALL CHECK_BLOWUP




  
   END DO

   IF(NumberStations>0)THEN

       CALL STATIONS

   ENDIF






     ! record wall time at the end

     if(myid.eq.0) tend = MPI_Wtime( )



     if(myid.eq.0) write(*,*) 'Simulation takes',tend-tbegin,'seconds'
     if(myid.eq.0) write(3,*) 'Simulation takes',tend-tbegin,'seconds'
     if (myid.eq.0) WRITE(*,*)'Normal Termination!'
     if (myid.eq.0) WRITE(3,*)'Normal Termination!'



     call MPI_FINALIZE ( ier )


END PROGRAM FUNWAVE_TVD


!-------------------------------------------------------------------------------------
! This part is not subroutines
!  DEFINITIONS OF VARIABLES
! 
!    Last Update: 02/18/2016 Fengyan Shi
!-------------------------------------------------------------------------------------
!
! Depth(): still water depth at element point
! DepthNode(): still water depth at node
! DepthX(): still water depth at x-interface
! DepthY(): still water depth at y-interface
! Eta():   surface elevation
! Eta0(): Eta at previous time level
!  for dry point, Eta() = MinDepth+Z()
! MASK(): 1 - wet
!         0 - dry
! MASK_STRUC(): 0 - permanent dry point
! MASK9: mask for itself and 8 elements around
! 
! U():  depth-averaged u or u at the reference level (u_alpha) at element
! V():  depth-averaged v or v at the reference level (v_alpha) at element
! HU(): (dep+eta)*u at element
! HV(): (dep+eta)*v at element
! P(): HU + dispersion at x-interface
! Q(): HV + dispersion at y-interface
! Fx(): F at x-interface
! Fy(): F at y-interface
! Gx(): G at x-interface
! Gy(): G at y-interface
! Ubar(:,:,:): Ubar
! Vbar(:,:,:): Vbar

! dispersion
! U1p(:,:): x-component of V1p
! V1p(:,:): y-component of V1p

! 
! EtaRxL(): Eta Left value at x-interface
! EtaRxR(): Eta Right value at x-interface
! EtaRyL(): Eta Left value at y-interface
! EtaRyR(): Eta Right value at y-interface
! HxL():   total depth  Left value at x-interface
! HxR():   total depth  Right value at x-interface
! HyL():   total depth  Left value at y-interface
! HyR():   total depth  Right value at y-interface

! HUxL(): HU Left value at x-interface
! HUxR(): HU Right value at x-interface
! HUyL(): HV Left value at y-interface
! HUyR(): HV Right value at y-interface

! PL(): HU + dispersion, Left value at x-interface
! PR(): HU + dispersion, Right value at x-interface
! QL(): HV + dispersion, Left value at y-interface
! QR(): HV + dispersion, Right value at y-interface

! FxL = HUxL*UxL + 1/2*g*(EtaRxL^2 + 2*EtaRxL*Depthx)
! FxR = HUxR*UxR + 1/2*g*(EtaRxR^2 + 2*EtaRxR*Depthx)
! FyL = HyL*UyL*VyL
! FyR = HyR*UyR*VyR

! GxL = HxL*UxL*VxL
! GxR = HxR*UxR*VxR
! GyL = HVyL*VyL + 1/2*g*(EtaRyL^2 + 2*EtaRyL*Depthy)
! GyR = HVyR*VyR + 1/2*g*(EtaRyR^2 + 2*EtaRyR*Depthy) 






