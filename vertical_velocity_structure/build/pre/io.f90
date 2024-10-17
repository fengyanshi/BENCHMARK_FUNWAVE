!------------------------------------------------------------------------------------
!
!      FILE io.F
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
!
!    OUTPUT is subroutine for screen, station, and field print-out
!
!    HISTORY:
!      12/06/2017  Michael-Angelo Y.-H. Lam
!      01/10/2011  Fengyan SHi
!-------------------------------------------------------------------------------------
SUBROUTINE OUTPUT
    USE GLOBAL


    IMPLICIT NONE

     SCREEN_COUNT=SCREEN_COUNT+DT

     IF(SCREEN_COUNT>=SCREEN_INTV)THEN
      SCREEN_COUNT=SCREEN_COUNT-SCREEN_INTV
      CALL STATISTICS
     ENDIF

! stations
      IF(NumberStations>0)THEN
      PLOT_COUNT_STATION=PLOT_COUNT_STATION+DT
      IF(PLOT_COUNT_STATION>=PLOT_INTV_STATION)THEN
       PLOT_COUNT_STATION=PLOT_COUNT_STATION-PLOT_INTV_STATION

       CALL STATIONS

      ENDIF
      ENDIF
! preview

      IF(TIME>=PLOT_START_TIME)THEN

	PLOT_COUNT=PLOT_COUNT+DT
      IF(PLOT_COUNT>=PLOT_INTV)THEN
       PLOT_COUNT=PLOT_COUNT-PLOT_INTV
       CALL PREVIEW
      ENDIF



   ENDIF ! end plot start time


END SUBROUTINE OUTPUT

!-------------------------------------------------------------------------------------
!
!    READ_INPUT is subroutine to read from input.txt
!
!  HISTORY:
!  01/10/2011  Fengyan SHi
!  12/23/2014  Young-Kwang Choi, added option for intel compiler
!  07/17/2019  Zhouteng Ye, added 1. get input from command line argument
!                                 2. input file for post-processor
!                       
!
!-------------------------------------------------------------------------------------

SUBROUTINE READ_INPUT
    USE GLOBAL
    USE INPUT_READ
!    USE TIME_SPECTRA_MODULE
![jychoi added this for intel compiler 14.12.23    


!jychoi 14.12.23]
    IMPLICIT NONE
    CHARACTER(LEN=80) FILE_NAME
    CHARACTER(LEN=80) MKFOLDER
    INTEGER::LINE
    INTEGER :: ierr
    INTEGER :: I_comp
    LOGICAL :: INPUT_PHASE = .FALSE.

    !>by Zhouteng Ye
    CHARACTER(LEN=80)::INPUT_NAME=''

![ykchoi
    CHARACTER(LEN=80)::FDIR=' '
!ykchoi]


    CALL MPI_COMM_SIZE (MPI_COMM_WORLD, nprocs, ier)   !ykchoi(04/May/2017)
    CALL MPI_COMM_RANK (MPI_COMM_WORLD, myid, ier)


      FDIR=TRIM(RESULT_FOLDER)
	OPEN(10000,FILE='time_dt.out',STATUS='UNKNOWN')

      OPEN(3,FILE='LOG.txt')   

! read everything from input.txt

      !>by Zhouteng Ye
      !> Get the argument from the command.
      !> If no input in command, file name is 'input.txt' (same as before)
      !> If the input comes with other name, read the correponding file
      CALL GETARG(1,INPUT_NAME) 
      if (INPUT_NAME .eq. '') Then
        FILE_NAME='input.txt'
      Else
        FILE_NAME=INPUT_NAME
      endif
      INPUT_FILE_NAME=FILE_NAME

! title
      CALL READ_STRING(TITLE,FILE_NAME,'TITLE',ierr)
      IF(ierr==1)THEN
        !write(*,*) 'No TITLE in ', FILE_NAME, 'use default'
        TITLE='---TEST RUN---'
      ENDIF

      if (myid.eq.0) WRITE(3,*)'-------------- LOG FILE -----------------'
      if (myid.eq.0) WRITE(3,*)TITLE
      if (myid.eq.0) WRITE(3,*)' --------------input start --------------'



      if (myid.eq.0) WRITE(3,*)'                                         '
      if (myid.eq.0) WRITE(3,*)'-------------- PARALLEL -----------------'    




! parallel info
      CALL READ_INTEGER(PX,FILE_NAME,'PX',ierr)
      IF(ierr == 1) THEN
        PX = 1
        if (myid.eq.0)write(*,*) 'No PX sepecified ', 'use PX=1'
        if (myid.eq.0)WRITE(3,'(A20,A20)')'No PX sepecified ', 'use PX=1'
      ENDIF
      CALL READ_INTEGER(PY,FILE_NAME,'PY',ierr)  
      IF(ierr == 1) THEN
        PY = 1
        if (myid.eq.0)write(*,*) 'No PY sepecified ', 'use PY=1'
        if (myid.eq.0)WRITE(3,'(A20,A20)')'No PY sepecified ', 'use PY=1'
      ENDIF
      if (myid.eq.0) WRITE(3,'(A7,I3,A7,I3)') 'PX   =',PX,'PY   =', PY




      if (myid.eq.0) WRITE(3,*)'                                         '
      if (myid.eq.0) WRITE(3,*)'-------------- GRID INFO -----------------'


! dimension
      CALL READ_INTEGER(Mglob,FILE_NAME,'Mglob',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'Mglob:', 'NOT DEFINED, STOP'
         WRITE(3,'(A40,A40)')'Mglob:', 'NOT DEFINED, STOP'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

      CALL READ_INTEGER(Nglob,FILE_NAME,'Nglob',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'Nglob:', 'NOT DEFINED, STOP'
         WRITE(3,'(A40,A40)')'Nglob:', 'NOT DEFINED, STOP'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF


      if (myid.eq.0) WRITE(3,'(A7,I8,A7,I8)') 'Mglob=',Mglob,'Nglob=', Nglob

! grid 


      CALL READ_FLOAT(DX,FILE_NAME,'DX',ierr)

      IF(ierr==1)THEN
         PRINT *,"Did you intend to use Spherical Coordinates?"

      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'DX:', 'NOT DEFINED, STOP'
         WRITE(3,'(A40,A40)')'DX:', 'NOT DEFINED, STOP'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF


      CALL READ_FLOAT(DY,FILE_NAME,'DY',ierr)

      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'DY:', 'NOT DEFINED, STOP'
         WRITE(3,'(A40,A40)')'DY:', 'NOT DEFINED, STOP'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF


      if (myid.eq.0) WRITE(3,'(A4,F12.2,A4,F12.2)')'DX=',DX,'DY=',DY



  ! end spherical

! depth 
      CALL READ_STRING(DEPTH_TYPE,FILE_NAME,'DEPTH_TYPE',ierr)

      IF(ierr==1)THEN
        DEPTH_TYPE = 'FLAT'

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'You dont specify DEPTH_TYPE, use FLAT.'
         WRITE(3,'(A40)')'You dont specify DEPTH_TYPE, use FLAT.'
      endif

      ENDIF


      if (myid.eq.0) WRITE(3,'(A12,A50)')'DEPTH_TYPE:', DEPTH_TYPE


      IF(DEPTH_TYPE(1:3)=='DAT')THEN
        CALL READ_STRING(DEPTH_FILE,FILE_NAME,'DEPTH_FILE',ierr)

      if (myid.eq.0) WRITE(3,'(A12,A50)')'DEPTH_FILE:', DEPTH_FILE

      ENDIF  ! end type=data

    IF(DEPTH_TYPE(1:3)=='FLA')THEN
      CALL READ_FLOAT(DEPTH_FLAT,FILE_NAME,'DEPTH_FLAT',ierr) 
      
      IF(ierr==1)THEN
        DEPTH_FLAT = 10.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'You dont specify DEPTH_FLAT, use 10 m.'
         WRITE(3,'(A40)')'You dont specify DEPTH_FLAT, use 10 m.'
      endif

      ENDIF      


      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'DEPTH_FLAT=', DEPTH_FLAT  

    ENDIF ! endif type=flat


      IF(DEPTH_TYPE(1:3)=='SLO')THEN
      CALL READ_FLOAT(DEPTH_FLAT,FILE_NAME,'DEPTH_FLAT',ierr) 

      IF(ierr==1)THEN
        DEPTH_FLAT = 10.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'You dont specify DEPTH_FLAT, use 10 m.'
         WRITE(3,'(A40)')'You dont specify DEPTH_FLAT, use 10 m.'
      endif

      ENDIF 

      CALL READ_FLOAT(SLP,FILE_NAME,'SLP',ierr) 

      IF(ierr==1)THEN
        SLP = 0.1_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'You dont specify SLP, use 0.1'
         WRITE(3,'(A40)')'You dont specify SLP, use 0.1'
      endif

      ENDIF 

      CALL READ_FLOAT(Xslp,FILE_NAME,'Xslp',ierr) 

      IF(ierr==1)THEN
        Xslp = 0.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'You dont specify Xslp, use 0.0'
         WRITE(3,'(A40)')'You dont specify Xslp, use 0.0'
      endif

      ENDIF 


      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'DEPTH_FLAT=', DEPTH_FLAT 
      if (myid.eq.0) WRITE(3,'(A5,F12.2)')'SLP=', SLP
      if (myid.eq.0) WRITE(3,'(A6,F12.2)')'Xslp=', Xslp  

      ENDIF  ! endif type=slope

! depth correction

      CALL READ_LOGICAL(BATHY_CORRECTION,FILE_NAME,'BATHY_CORRECTION',ierr)
      IF(ierr == 1)THEN
       BATHY_CORRECTION = .FALSE. 
      ENDIF
      IF(BATHY_CORRECTION)THEN

      if (myid.eq.0)then
       WRITE(3,'(A40)')'Bathymetry is corrected !'
       WRITE(*,'(A40)')'Bathymetry is corrected !'
      endif

      ENDIF

! time


      if (myid.eq.0) WRITE(3,*)'                                         '
      if (myid.eq.0) WRITE(3,*)'-------------- TIME INFO -----------------'


      CALL READ_FLOAT(TOTAL_TIME,FILE_NAME,'TOTAL_TIME',ierr)

      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'TOTAL_TIME:', 'NOT FOUND, STOP'
         WRITE(3,'(A40,A40)')'TOTAL_TIME:', 'NOT FOUND, STOP'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

      CALL READ_FLOAT(PLOT_START_TIME,FILE_NAME,'PLOT_START_TIME',ierr)

      IF(ierr==1)THEN
        PLOT_START_TIME = 0.0

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'PLOT_START_TIME Default:  0.0 s'
         WRITE(3,'(A40)')'PLOT_START_TIME Default:  0.0 s'
      endif

       ENDIF

      CALL READ_FLOAT(PLOT_INTV,FILE_NAME,'PLOT_INTV',ierr)

      IF(ierr==1)THEN
        PLOT_INTV = 1.0

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'PLOT_INTV Default:  1.0 s'
         WRITE(3,'(A40)')'PLOT_INTV Default:  1.0 s'
      endif

       ENDIF

      CALL READ_FLOAT(PLOT_INTV_STATION,FILE_NAME,'PLOT_INTV_STATION',ierr)

      IF(ierr==1)THEN
        PLOT_INTV_STATION = 1.0

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'PLOT_INTV_STATION Default:  1.0 s'
         WRITE(3,'(A40)')'PLOT_INTV_STATION Default:  1.0 s'
      endif

       ENDIF

      CALL READ_INTEGER(StationOutputBuffer,FILE_NAME,'StationOutputBuffer',ierr)
      IF(ierr==1)THEN
        StationOutputBuffer = 1000

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'StationOutputBuffer not specified, use default:1000'
         WRITE(3,'(A80)')'StationOutputBuffer not specified, use default:1000'
      endif

      ENDIF

      CALL READ_FLOAT(SCREEN_INTV,FILE_NAME,'SCREEN_INTV',ierr)
      IF(ierr==1)THEN
        SCREEN_INTV = 1.0

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'SCREEN_INTV Default:  1.0 s'
         WRITE(3,'(A40)')'SCREEN_INTV Default:  1.0 s'
      endif

       ENDIF



      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'TOTAL_TIME=', TOTAL_TIME
      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'PLOT_INTV= ', PLOT_INTV
      if (myid.eq.0) WRITE(3,'(A13,F12.2)')'SCREEN_INTV=', SCREEN_INTV



      if (myid.eq.0) WRITE(3,*)'                                         '
      if (myid.eq.0) WRITE(3,*)'-------------- HOT START -----------------'


! initial uvz
      CALL READ_LOGICAL(INI_UVZ,FILE_NAME,'INI_UVZ',ierr)
      IF(ierr==1)THEN
        INI_UVZ = .FALSE.
      ENDIF
    IF(INI_UVZ)THEN

        CALL READ_LOGICAL(BED_DEFORMATION,FILE_NAME,'BED_DEFORMATION',ierr)
        IF(ierr==1)THEN
           BED_DEFORMATION = .FALSE.
        ENDIF 

        IF(BED_DEFORMATION)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'Bathymetry is adjusted based on Bed deformation.'
         WRITE(3,'(A50)')'Bathymetry is adjusted based on Bed deformation.'
      endif

        ENDIF ! end bed deformation

        CALL READ_STRING(ETA_FILE,FILE_NAME,'ETA_FILE',ierr)

      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'You use INI_UVZ, ETA_FILE NOT FOUND, STOP'
         WRITE(3,'(A50)')'You use INI_UVZ, ETA_FILE NOT FOUND, STOP'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

        CALL READ_STRING(U_FILE,FILE_NAME,'U_FILE',ierr)

        IF(ierr==1)THEN
          NO_UV_FILE = .TRUE.

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'U_FILE Default:  ZERO.'
         WRITE(3,'(A40)')'U_FILE Default:  ZERO.'
      endif

        ELSE
          NO_UV_FILE = .FALSE.
        ENDIF

        CALL READ_STRING(V_FILE,FILE_NAME,'V_FILE',ierr)

        IF(ierr==1)THEN
          NO_UV_FILE = .TRUE.

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'V_FILE Default:  ZERO.'
         WRITE(3,'(A40)')'V_FILE Default:  ZERO.'
      endif

        ELSE
          NO_UV_FILE = .FALSE.
          IF(NO_UV_FILE) NO_UV_FILE = .TRUE.    ! in case you dont have u file
        ENDIF

        CALL READ_STRING(MASK_FILE,FILE_NAME,'MASK_FILE',ierr)

        IF(ierr==1)THEN
          NO_MASK_FILE = .TRUE.

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'MASK_FILE NOT FOUND, USE NO_MASK option'
         WRITE(3,'(A40)')'MASK_FILE NOT FOUND, USE NO_MASK option'
      endif


        ELSE
          NO_MASK_FILE = .FALSE.
        ENDIF

	![ykchoi(14.12.24.)
	  CALL READ_FLOAT(HotStartTime,FILE_NAME,'HotStartTime',ierr)

      IF(ierr==1)THEN
        HotStartTime = 0.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'HotStartTime Default:  0.0 s'
         WRITE(3,'(A40)')'HotStartTime Default:  0.0 s'
      endif

       ENDIF

	!ykchoi(14.12.24.)]
        CALL READ_INTEGER(icount,FILE_NAME,'OutputStartNumber',ierr)
      IF(ierr==1)THEN
        icount = 1

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'OutputStartNumber Default:  1'
         WRITE(3,'(A40)')'OutputStartNumber Default:  1'
      endif

       ENDIF

        icount = icount-1

       ENDIF  ! end initial hot start file


      if (myid.eq.0) WRITE(3,*)'                                         '
      if (myid.eq.0) WRITE(3,*)'-------------- STILL WATER INFO -----------------'


! add water level 03/29/2016

        CALL READ_FLOAT(WaterLevel,FILE_NAME,'WaterLevel',ierr)
        IF(ierr==1)THEN
          WaterLevel = 0.0
        ENDIF


      if (myid.eq.0) WRITE(3,'(A20,F12.5)')'WaterLevel = ', WaterLevel




      if (myid.eq.0) WRITE(3,*)'                                         '
      if (myid.eq.0) WRITE(3,*)'-------------- WAVEMAKER -----------------'


! wavemaker
      CALL READ_STRING(WaveMaker,FILE_NAME,'WAVEMAKER',ierr)
      IF(ierr==1)THEN
        WaveMaker = 'nothing'

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'No WaveMaker'
         WRITE(3,'(A40)')'No WaveMaker'
      endif

       ENDIF


      if (myid.eq.0) WRITE(3,'(A11,A50)')'WAVEMAKER:', WAVEMAKER


        IF(WaveMaker(1:7)=='LEF_SOL')THEN
          CALL READ_FLOAT(AMP_SOLI,FILE_NAME,'AMP',ierr)

      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'AMP_SOLI NOT FOUND, specify AMP in input.txt'
         WRITE(3,'(A60)')'AMP_SOLI NOT FOUND, specify AMP in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(DEP_SOLI,FILE_NAME,'DEP',ierr)

      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'DEP_SOLI NOT FOUND, specify DEP in input.txt'
         WRITE(3,'(A60)')'DEP_SOLI NOT FOUND, specify DEP in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF


          CALL READ_FLOAT(LAG_SOLI,FILE_NAME,'LAGTIME',ierr)

      IF(ierr==1)THEN
       LAG_SOLI = 0.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'LAGTIME Default:  0.0'
         WRITE(3,'(A40)')'LAGTIME Default:  0.0'
      endif

       ENDIF


      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'AMP_SOLI=', AMP_SOLI
      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'DEP_SOLI=', DEP_SOLI
      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'LAG_SOLI=', LAG_SOLI

        ENDIF

        IF(WaveMaker(1:7)=='WK_TIME')THEN
        CALL READ_INTEGER(NumWaveComp,FILE_NAME,'NumWaveComp',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'NumWaveComp NOT FOUND, specify NumWaveComp in input.txt'
         WRITE(3,'(A80)')'NumWaveComp NOT FOUND, specify NumWaveComp in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

        CALL READ_FLOAT(PeakPeriod,FILE_NAME,'PeakPeriod',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'PeakPeriod NOT FOUND, specify PeakPeriod in input.txt'
         WRITE(3,'(A80)')'PeakPeriod NOT FOUND, specify PeakPeriod in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

        CALL READ_STRING(WaveCompFile,FILE_NAME,'WaveCompFile',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'WaveCompFile NOT FOUND, specify WaveCompFile in input.txt'
         WRITE(3,'(A80)')'WaveCompFile NOT FOUND, specify WaveCompFile in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(Xc_WK,FILE_NAME,'Xc_WK',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'Xc_WK NOT FOUND, specify Xc_WK in input.txt'
         WRITE(3,'(A80)')'Xc_WK NOT FOUND, specify Xc_WK in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(Yc_WK,FILE_NAME,'Yc_WK',ierr)
      IF(ierr==1)THEN
       Yc_WK = ZERO

      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'Yc_WK defalt: 0.0'
         WRITE(3,'(A50)')'Yc_WK defalt: 0.0'
      endif

       ENDIF

          CALL READ_FLOAT(DEP_WK,FILE_NAME,'DEP_WK',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'DEP_WK NOT FOUND, specify DEP_WK in input.txt'
         WRITE(3,'(A80)')'DEP_WK NOT FOUND, specify DEP_WK in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(Time_ramp,FILE_NAME,'Time_ramp',ierr)
      IF(ierr==1)THEN
        Time_ramp = 0.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Time_ramp Default:  0.0'
         WRITE(3,'(A40)')'Time_ramp Default:  0.0'
      endif

       ENDIF

          CALL READ_FLOAT(Delta_WK,FILE_NAME,'Delta_WK',ierr)
      IF(ierr==1)THEN
        Delta_WK = 0.5_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Delta_WK Default:  0.5'
         WRITE(3,'(A40)')'Delta_WK Default:  0.5'
      endif

       ENDIF


          CALL READ_FLOAT(Ywidth_WK,FILE_NAME,'Ywidth_WK',ierr)
      IF(ierr==1)THEN
        Ywidth_WK = LARGE

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Ywidth_WK Default:  LARGE'
         WRITE(3,'(A40)')'Ywidth_WK Default:  LARGE'
      endif

       ENDIF


      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'Xc_WK   =', Xc_WK
      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'Yc_WK   =', Yc_WK
      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'DEP_WK  =', DEP_WK
      if (myid.eq.0) WRITE(3,'(A11,F12.2)')'Time_ramp=', Time_ramp
      if (myid.eq.0) WRITE(3,'(A11,F12.2)')'Delta_WK=', Delta_WK
      if (myid.eq.0) WRITE(3,'(A11,F12.2)')'Ywidth_WK=', Ywidth_WK


        ENDIF  ! end WK_TIME

        IF(WaveMaker(1:7)=='INI_SOL')THEN

          CALL READ_LOGICAL(SolitaryPositiveDirection,FILE_NAME,  &
               'SolitaryPositiveDirection',ierr)

        IF(ierr==1)THEN
          SolitaryPositiveDirection = .TRUE.
        ENDIF

        IF(SolitaryPositiveDirection) THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'Solitary wave propagate in + X direction'
         WRITE(3,'(A60)')'Solitary wave propagate in + X direction'
      endif

        ELSE

      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'Solitary wave propagate in - X direction'
         WRITE(3,'(A60)')'Solitary wave propagate in - X direction'
      endif

        ENDIF

          CALL READ_FLOAT(AMP_SOLI,FILE_NAME,'AMP',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'AMP_SOLI NOT FOUND, specify AMP in input.txt'
         WRITE(3,'(A60)')'AMP_SOLI NOT FOUND, specify AMP in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(DEP_SOLI,FILE_NAME,'DEP',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'DEP_SOLI NOT FOUND, specify DEP in input.txt'
         WRITE(3,'(A60)')'DEP_SOLI NOT FOUND, specify DEP in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(XWAVEMAKER,FILE_NAME,'XWAVEMAKER',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'XWAVEMAKER NOT FOUND, specify XWAVEMAKER in input.txt'
         WRITE(3,'(A80)')'XWAVEMAKER NOT FOUND, specify XWAVEMAKER in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF


      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'AMP_SOLI=', AMP_SOLI
      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'DEP_SOLI=', DEP_SOLI

        ENDIF  ! end initial solitary

        IF(WaveMaker(1:6)=='N_WAVE')THEN
          CALL READ_FLOAT(x1_Nwave,FILE_NAME,'x1_Nwave',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'x1_Nwave NOT FOUND, specify x1_Nwave in input.txt'
         WRITE(3,'(A80)')'x1_Nwave NOT FOUND, specify x1_Nwave in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(x2_Nwave,FILE_NAME,'x2_Nwave',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'x2_Nwave NOT FOUND, specify x2_Nwave in input.txt'
         WRITE(3,'(A80)')'x2_Nwave NOT FOUND, specify x2_Nwave in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(a0_Nwave,FILE_NAME,'a0_Nwave',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'a0_Nwave NOT FOUND, specify a0_Nwave in input.txt'
         WRITE(3,'(A80)')'a0_Nwave NOT FOUND, specify a0_Nwave in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(gamma_Nwave,FILE_NAME,'gamma_Nwave',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'gamma_Nwave NOT FOUND, specify gamma_Nwave in input.txt'
         WRITE(3,'(A80)')'gamma_Nwave NOT FOUND, specify gamma_Nwave in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(dep_Nwave,FILE_NAME,'dep_Nwave',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'dep_Nwave NOT FOUND, specify dep_Nwave in input.txt'
         WRITE(3,'(A80)')'dep_Nwave NOT FOUND, specify dep_Nwave in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF


      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'x1_Nwave=', x1_Nwave
      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'x2_Nwave=', x2_Nwave
      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'a0_Nwave=', a0_Nwave
      if (myid.eq.0) WRITE(3,'(A13,F12.2)')'gamma_Nwave=', gamma_Nwave
      if (myid.eq.0) WRITE(3,'(A11,F12.2)')'dep_Nwave=', dep_Nwave

        ENDIF  ! end N_wave

        IF(WaveMaker(1:7)=='INI_REC')THEN
          CALL READ_FLOAT(AMP_SOLI,FILE_NAME,'AMP',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'AMP NOT FOUND, specify AMP in input.txt'
         WRITE(3,'(A50)')'AMP NOT FOUND, specify AMP in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(Xc,FILE_NAME,'Xc',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Xc NOT FOUND, specify Xc in input.txt'
         WRITE(3,'(A40)')'Xc NOT FOUND, specify Xc in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(Yc,FILE_NAME,'Yc',ierr)
      IF(ierr==1)THEN
       Yc = ZERO

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Yc NOT FOUND, specify Yc in input.txt'
         WRITE(3,'(A40)')'Yc NOT FOUND, specify Yc in input.txt'
      endif

      ENDIF


          CALL READ_FLOAT(WID,FILE_NAME,'WID',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'WID NOT FOUND, specify WID in input.txt'
         WRITE(3,'(A50)')'WID NOT FOUND, specify WID in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF


      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'AMP     =', AMP_SOLI
      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'Xc      =', Xc
      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'Yc      =', Yc
      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'WID     =', WID

        ENDIF ! endif rectangular hump

        IF(WaveMaker(1:7)=='INI_GAU'.OR.&
           WaveMaker(1:7)=='INI_DIP')THEN

          CALL READ_FLOAT(AMP_SOLI,FILE_NAME,'AMP',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'AMP NOT FOUND, specify AMP in input.txt'
         WRITE(3,'(A50)')'AMP NOT FOUND, specify AMP in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(Xc,FILE_NAME,'Xc',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'Xc NOT FOUND, specify Xc in input.txt'
         WRITE(3,'(A50)')'Xc NOT FOUND, specify Xc in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(Yc,FILE_NAME,'Yc',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'Yc NOT FOUND, specify Yc in input.txt'
         WRITE(3,'(A50)')'Yc NOT FOUND, specify Yc in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(WID,FILE_NAME,'WID',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'WID NOT FOUND, specify WID in input.txt'
         WRITE(3,'(A50)')'WID NOT FOUND, specify WID in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF


      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'AMP     =', AMP_SOLI
      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'Xc      =', Xc
      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'Yc      =', Yc
      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'WID(gamma)=', WID

        ENDIF ! endif gaussian hump

        IF(WaveMaker(1:6)=='WK_REG')THEN
          CALL READ_FLOAT(Xc_WK,FILE_NAME,'Xc_WK',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'Xc_WK NOT FOUND, specify Xc_WK in input.txt'
         WRITE(3,'(A60)')'Xc_WK NOT FOUND, specify Xc_WK in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(Yc_WK,FILE_NAME,'Yc_WK',ierr)
      IF(ierr==1)THEN
       Yc_WK = ZERO

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Yc_WK defalt: 0.0'
         WRITE(3,'(A40)')'Yc_WK defalt: 0.0'
      endif

       ENDIF


          CALL READ_FLOAT(Tperiod,FILE_NAME,'Tperiod',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'Tperiod NOT FOUND, specify Tperiod in input.txt'
         WRITE(3,'(A60)')'Tperiod NOT FOUND, specify Tperiod in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(AMP_WK,FILE_NAME,'AMP_WK',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'AMP_WK NOT FOUND, specify AMP_WK in input.txt'
         WRITE(3,'(A60)')'AMP_WK NOT FOUND, specify AMP_WK in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(DEP_WK,FILE_NAME,'DEP_WK',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'DEP_WK NOT FOUND, specify DEP_WK in input.txt'
         WRITE(3,'(A60)')'DEP_WK NOT FOUND, specify DEP_WK in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(Theta_WK,FILE_NAME,'Theta_WK',ierr)
      IF(ierr==1)THEN
        Theta_WK = 0.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Theta_WK Default:  0.0'
         WRITE(3,'(A40)')'Theta_WK Default:  0.0'
      endif

       ENDIF

          CALL READ_FLOAT(Time_ramp,FILE_NAME,'Time_ramp',ierr)
      IF(ierr==1)THEN
        Time_ramp = 0.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Time_ramp Default:  0.0'
         WRITE(3,'(A40)')'Time_ramp Default:  0.0'
      endif

       ENDIF

          CALL READ_FLOAT(Delta_WK,FILE_NAME,'Delta_WK',ierr)
      IF(ierr==1)THEN
        Delta_WK = 0.5_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Delta_WK Default:  0.5'
         WRITE(3,'(A40)')'Delta_WK Default:  0.5'
      endif

       ENDIF

          CALL READ_FLOAT(Ywidth_WK,FILE_NAME,'Ywidth_WK',ierr)
      IF(ierr==1)THEN
        Ywidth_WK = LARGE

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Ywidth_WK Default:  LARGE'
         WRITE(3,'(A40)')'Ywidth_WK Default:  LARGE'
      endif

       ENDIF


      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'Xc_WK   =', Xc_WK
      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'Yc_WK   =', Yc_WK
      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'Tperiod =', Tperiod
      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'AMP_WK  =', AMP_WK
      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'DEP_WK  =', DEP_WK
      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'Theta_WK=', Theta_WK
      if (myid.eq.0) WRITE(3,'(A11,F12.2)')'Time_ramp=', Time_ramp
      if (myid.eq.0) WRITE(3,'(A11,F12.2)')'Delta_WK=', Delta_WK
      if (myid.eq.0) WRITE(3,'(A11,F12.2)')'Ywidth_WK=', Ywidth_WK

        ENDIF  ! endif WK_REG

        IF(WaveMaker(1:6)=='WK_IRR'.OR.WaveMaker(1:6)=='TMA_1D'  &
           .OR.WaveMaker(1:6)=='JON_1D'.OR.WaveMaker(1:6)=='JON_2D')THEN

          CALL READ_FLOAT(Xc_WK,FILE_NAME,'Xc_WK',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'Xc_WK NOT FOUND, specify Xc_WK in input.txt'
         WRITE(3,'(A60)')'Xc_WK NOT FOUND, specify Xc_WK in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(Yc_WK,FILE_NAME,'Yc_WK',ierr)
      IF(ierr==1)THEN
       Yc_WK = ZERO

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Yc_WK defalt: 0.0'
         WRITE(3,'(A40)')'Yc_WK defalt: 0.0'
      endif

      ENDIF

          CALL READ_FLOAT(DEP_WK,FILE_NAME,'DEP_WK',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'DEP_WK NOT FOUND, specify DEP_WK in input.txt'
         WRITE(3,'(A60)')'DEP_WK NOT FOUND, specify DEP_WK in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(Time_ramp,FILE_NAME,'Time_ramp',ierr)
      IF(ierr==1)THEN
        Time_ramp = 0.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Time_ramp Default:  0.0'
         WRITE(3,'(A40)')'Time_ramp Default:  0.0'
      endif

       ENDIF

          CALL READ_FLOAT(Delta_WK,FILE_NAME,'Delta_WK',ierr)
      IF(ierr==1)THEN
        Delta_WK = 0.5_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Delta_WK Default:  0.5'
         WRITE(3,'(A40)')'Delta_WK Default:  0.5'
      endif

       ENDIF

          CALL READ_FLOAT(FreqPeak,FILE_NAME,'FreqPeak',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'FreqPeak NOT FOUND, specify FreqPeak in input.txt'
         WRITE(3,'(A80)')'FreqPeak NOT FOUND, specify FreqPeak in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(FreqMin,FILE_NAME,'FreqMin',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'FreqMin NOT FOUND, specify FreqMin in input.txt'
         WRITE(3,'(A80)')'FreqMin NOT FOUND, specify FreqMin in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(FreqMax,FILE_NAME,'FreqMax',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'FreqMax NOT FOUND, specify FreqMax in input.txt'
         WRITE(3,'(A80)')'FreqMax NOT FOUND, specify FreqMax in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(Hmo,FILE_NAME,'Hmo',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'Hmo NOT FOUND, specify Hmo in input.txt'
         WRITE(3,'(A50)')'Hmo NOT FOUND, specify Hmo in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF
          CALL READ_FLOAT(GammaTMA,FILE_NAME,'GammaTMA',ierr)
          IF(ierr==1)THEN
             GammaTMA = 3.3

      if (myid.eq.0) WRITE(3,'(A14,A50)')'GammaTMA', 'NOT DEFINED, USE 3.3'

          ENDIF

          CALL READ_INTEGER(Nfreq,FILE_NAME,'Nfreq',ierr)
          IF(ierr==1)THEN
             Nfreq = 45

      if (myid.eq.0) WRITE(3,'(A25,A50)')'Nfreq for TMA or JON', 'NOT DEFINED, USE 45'

          ENDIF

         IF(WaveMaker(1:6)=='TMA_1D'  &
           .OR.WaveMaker(1:6)=='JON_1D')THEN
          Ntheta = 1
         ELSE
          CALL READ_INTEGER(Ntheta,FILE_NAME,'Ntheta',ierr)
          IF(ierr==1)THEN
             Ntheta = 24

      if (myid.eq.0) WRITE(3,'(A25,A50)')'Ntheta for TMA or JON', 'NOT DEFINED, USE 24'

          ENDIF

         ENDIF ! 1D or not

!         IF(WaveMaker(1:6)=='TMA_1D'  &
!           .OR.WaveMaker(1:6)=='JON_1D')THEN
!          ThetaPeak = 0.0_SP
!         ELSE

          CALL READ_FLOAT(ThetaPeak,FILE_NAME,'ThetaPeak',ierr)
          IF(ierr==1)THEN
             ThetaPeak = 0.0_SP

      if (myid.eq.0) WRITE(3,'(A25,A50)')'ThetaPeak', 'NOT DEFINED, USE 0.0'

          ENDIF ! end ierr

!         ENDIF ! 1D or not

         IF(WaveMaker(1:6)=='TMA_1D'  &
           .OR.WaveMaker(1:6)=='JON_1D')THEN
          ! do nothing
         ELSE
          CALL READ_FLOAT(Sigma_Theta,FILE_NAME,'Sigma_Theta',ierr)
          IF(ierr==1)THEN
             Sigma_Theta = 10.0_SP

      if (myid.eq.0) WRITE(3,'(A25,A50)')'Sigma_Theta', 'NOT DEFINED, USE 10.0'

          ENDIF ! end ierr

         ENDIF ! 1D or not

          CALL READ_FLOAT(Ywidth_WK,FILE_NAME,'Ywidth_WK',ierr)
      IF(ierr==1)THEN
        Ywidth_WK = LARGE

      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'Ywidth_WK Default:  LARGE'
         WRITE(3,'(A50)')'Ywidth_WK Default:  LARGE'
      endif

       ENDIF



      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'Xc_WK   =  ', Xc_WK
      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'Yc_WK   =', Yc_WK
      if (myid.eq.0) WRITE(3,'(A11,F12.2)')'Ywidth_WK=', Ywidth_WK
      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'DEP_WK  =  ', DEP_WK
      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'Time_ramp= ', Time_ramp
      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'Delta_WK=  ', Delta_WK
      if (myid.eq.0) WRITE(3,'(A12,I12)')'Nfreq=  ', Nfreq
      if (myid.eq.0) WRITE(3,'(A12,I12)')'Ntheta=  ', Ntheta
      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'FreqPeak=  ', FreqPeak
      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'FreqMin =  ', FreqMin
      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'FreqMax =  ', FreqMax
      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'Hmo     =  ', Hmo
      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'GammaTMA=  ', GammaTMA
      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'ThetaPeak= ', ThetaPeak
      if (myid.eq.0) WRITE(3,'(A13,F12.2)')'Sigma_Theta=', Sigma_Theta

        ENDIF ! endif wk_irr 

       CALL READ_LOGICAL(ETA_LIMITER,FILE_NAME,'ETA_LIMITER',ierr)
       IF(ETA_LIMITER)THEN
          CALL READ_FLOAT(CrestLimit,FILE_NAME,'CrestLimit',ierr)
          CALL READ_FLOAT(TroughLimit,FILE_NAME,'TroughLimit',ierr)
       ENDIF
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!START!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Oct, 2021
!Salatin, R., Chen, Q., Bak, A. S., Shi, F., & Brandt, S. R. (2021). Effects of
!wave coherence on longshore variability of nearshore wave processes. Journal
!of Geophysical Research: Oceans,126, e2021JC017641.
!https://doi.org/10.1029/2021JC017641

    IF(WaveMaker(1:13)=='WK_NEW_DATA2D'.OR.WaveMaker(1:9)=='WK_DATA2D')THEN
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          CALL READ_FLOAT(Xc_WK,FILE_NAME,'Xc_WK',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'Xc_WK NOT FOUND, specify Xc_WK in input.txt'
         WRITE(3,'(A80)')'Xc_WK NOT FOUND, specify Xc_WK in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(Yc_WK,FILE_NAME,'Yc_WK',ierr)
      IF(ierr==1)THEN
         Yc_WK = ZERO

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'Yc_WK NOT FOUND, specify Yc_WK in input.txt'
         WRITE(3,'(A80)')'Yc_WK NOT FOUND, specify Yc_WK in input.txt'
      endif

      ENDIF

          CALL READ_FLOAT(DEP_WK,FILE_NAME,'DEP_WK',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'DEP_WK NOT FOUND, specify DEP_WK in input.txt'
         WRITE(3,'(A80)')'DEP_WK NOT FOUND, specify DEP_WK in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(Time_ramp,FILE_NAME,'Time_ramp',ierr)
      IF(ierr==1)THEN
        Time_ramp = 0.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'Time_ramp Default:  0.0'
         WRITE(3,'(A50)')'Time_ramp Default:  0.0'
      endif

       ENDIF

          CALL READ_FLOAT(Delta_WK,FILE_NAME,'Delta_WK',ierr)
      IF(ierr==1)THEN
        Delta_WK = 0.5_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'Delta_WK Default:  0.5'
         WRITE(3,'(A50)')'Delta_WK Default:  0.5'
      endif

       ENDIF

          CALL READ_STRING(WaveCompFile,FILE_NAME,'WaveCompFile',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'WaveCompFile NOT FOUND, specify WaveCompFile in input.txt'
         WRITE(3,'(A80)')'WaveCompFile NOT FOUND, specify WaveCompFile in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

! data reading will be in wavemaker.F

          CALL READ_FLOAT(Ywidth_WK,FILE_NAME,'Ywidth_WK',ierr)
      IF(ierr==1)THEN
        Ywidth_WK = LARGE

      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'Ywidth_WK Default:  LARGE'
         WRITE(3,'(A50)')'Ywidth_WK Default:  LARGE'
      endif

       ENDIF


      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'Xc_WK   =  ', Xc_WK
      if (myid.eq.0) WRITE(3,'(A10,F12.2)')'Yc_WK   =', Yc_WK
      if (myid.eq.0) WRITE(3,'(A11,F12.2)')'Ywidth_WK=', Ywidth_WK
      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'DEP_WK  =  ', DEP_WK
      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'Time_ramp= ', Time_ramp
      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'Delta_WK=  ', Delta_WK

        ENDIF ! end ierr

       CALL READ_LOGICAL(EqualEnergy,FILE_NAME,'EqualEnergy',ierr)
        IF(ierr==1)THEN
             EqualEnergy = .FALSE.

      if (myid.eq.0) WRITE(3,'(A50)') 'EqualEnergy NOT USED in frequency domain'

          ENDIF 

! absorbing generating wavemaker, more tests needed 
     IF(WaveMaker(1:3)=='ABS'.OR.WaveMaker(1:11)=='LEFT_BC_IRR') THEN

       CALL READ_STRING(WAVE_DATA_TYPE,FILE_NAME,'WAVE_DATA_TYPE',ierr)

       CALL READ_FLOAT(DEP_Ser,FILE_NAME,'DepthWaveMaker',ierr)

        IF(ierr==1)THEN
       CALL READ_FLOAT(DEP_Ser,FILE_NAME,'DEP_WK',ierr)


      if (myid.eq.0) WRITE(3,'(A50)') 'DepthWaveMaker not defined, read DEP_WK'

      IF(ierr==1)THEN ! if can not find again

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'you should specify either DepthWaveMaker or DEP_WK,STOP'
         WRITE(3,'(A80)')'you should specify either DepthWaveMaker or DEP_WK,STOP'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          ENDIF ! end read DepthWaveMaker


      IF(WaveMaker(1:3)=='ABS')THEN
       CALL READ_FLOAT(WidthWaveMaker,FILE_NAME,'WidthWaveMaker',ierr) 
!  get R and A for sponge_wavemaker calculation  
        CALL READ_FLOAT(R_sponge_wavemaker,FILE_NAME,'R_sponge_wavemaker',ierr)
        CALL READ_FLOAT(A_sponge_wavemaker,FILE_NAME,'A_sponge_wavemaker',ierr)
      ENDIF

      IF(WAVE_DATA_TYPE(1:4)=='DATA')THEN
       CALL READ_STRING(WaveCompFile,FILE_NAME,'WaveCompFile',ierr)

      OPEN(1,FILE=TRIM(WaveCompFile))
       READ(1,*)NumFreq,NumDir
       ALLOCATE (Amp_Ser(NumFreq,NumDir),  &
          Per_Ser(NumFreq),Theta_Ser(NumDir),Phase_LEFT(NumFreq,NumDir))
       READ(1,*)PeakPeriod  ! useless for this application but should keep for consistency
       DO J=1,NumFreq
          READ(1,*)Per_Ser(J)  ! read in as frequency
!print*,J,Per_Ser(J)
       ENDDO
       DO I=1,NumDir
          READ(1,*)Theta_Ser(I)
       ENDDO
       DO I=1,NumDir
         READ(1,*)(Amp_Ser(J,I),J=1,NumFreq)
!print*,I,Amp_Ser(J,I)
       ENDDO
       DO I=1,NumDir
         READ(1,*,END=991)(Phase_LEFT(J,I),J=1,NumFreq)
       ENDDO
     CLOSE(1)
881  INPUT_PHASE = .TRUE.
991  CONTINUE

    DO J=1,NumFreq
     DO I=1,NumDir

     IF(INPUT_PHASE)THEN
         Phase_LEFT(J,I)=Phase_LEFT(J,I)*3.1415926/180.0_SP
     ELSE

          Phase_LEFT(J,I)=rand(0)*2.0_SP*3.1415926

     ENDIF

    ENDDO
   ENDDO

! to make consistent with cm and sm approach we use phase_ser which is 
! only is random with frequency
       ALLOCATE(Phase_Ser(NumFreq))

       ALLOCATE(Segma_Ser(NumFreq),Wave_Number_Ser(NumFreq) )
       DO J=1,NumFreq
          Phase_Ser(J)=Phase_LEFT(J,1)
          IF(Per_Ser(J).EQ.ZERO)THEN
          WRITE(*,*) 'wave frequency is zero, STOP'
          STOP
         ELSE
          Per_Ser(J)=1.0_SP/Per_Ser(J)
         ENDIF
       ENDDO
       DO I=1,NumDir
         Theta_Ser(I)=Theta_Ser(I)*DEG2RAD
       ENDDO  


      if (myid.eq.0) WRITE(3,'(A40)')'absorbing generating wave maker'
      if (myid.eq.0) WRITE(3,'(A40)')'use DATA'

 
     ELSE ! use TMA or JON

          CALL READ_FLOAT(FreqPeak,FILE_NAME,'FreqPeak',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'FreqPeak NOT FOUND, specify FreqPeak in input.txt'
         WRITE(3,'(A80)')'FreqPeak NOT FOUND, specify FreqPeak in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(FreqMin,FILE_NAME,'FreqMin',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'FreqMin NOT FOUND, specify FreqMin in input.txt'
         WRITE(3,'(A80)')'FreqMin NOT FOUND, specify FreqMin in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(FreqMax,FILE_NAME,'FreqMax',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'FreqMax NOT FOUND, specify FreqMax in input.txt'
         WRITE(3,'(A80)')'FreqMax NOT FOUND, specify FreqMax in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF

          CALL READ_FLOAT(Hmo,FILE_NAME,'Hmo',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'Hmo NOT FOUND, specify Hmo in input.txt'
         WRITE(3,'(A50)')'Hmo NOT FOUND, specify Hmo in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ENDIF
          CALL READ_FLOAT(GammaTMA,FILE_NAME,'GammaTMA',ierr)
          IF(ierr==1)THEN
             GammaTMA = 3.3

      if (myid.eq.0) WRITE(3,'(A14,A50)')'GammaTMA', 'NOT DEFINED, USE 3.3'

          ENDIF

          CALL READ_INTEGER(Nfreq,FILE_NAME,'Nfreq',ierr)
          IF(ierr==1)THEN
             Nfreq = 45

      if (myid.eq.0) WRITE(3,'(A25,A50)')'Nfreq for TMA or JON', 'NOT DEFINED, USE 45'

          ENDIF

         IF(WAVE_DATA_TYPE(1:6)=='TMA_1D'  &
           .OR.WAVE_DATA_TYPE(1:6)=='JON_1D')THEN
          Ntheta = 1
         ELSE
          CALL READ_INTEGER(Ntheta,FILE_NAME,'Ntheta',ierr)
          IF(ierr==1)THEN
             Ntheta = 24

      if (myid.eq.0) WRITE(3,'(A25,A50)')'Ntheta for TMA or JON', 'NOT DEFINED, USE 24'

          ENDIF

         ENDIF ! 1D or not

!         IF(WAVE_DATA_TYPE(1:6)=='TMA_1D'  &
!           .OR.WAVE_DATA_TYPE(1:6)=='JON_1D')THEN
!          ThetaPeak = 0.0_SP
!         ELSE

          CALL READ_FLOAT(ThetaPeak,FILE_NAME,'ThetaPeak',ierr)
          IF(ierr==1)THEN
             ThetaPeak = 0.0_SP

      if (myid.eq.0) WRITE(3,'(A25,A50)')'ThetaPeak', 'NOT DEFINED, USE 0.0'

          ENDIF ! end ierr

!         ENDIF ! 1D or not

         IF(WAVE_DATA_TYPE(1:6)=='TMA_1D'  &
           .OR.WAVE_DATA_TYPE(1:6)=='JON_1D')THEN
          ! do nothing
         ELSE
          CALL READ_FLOAT(Sigma_Theta,FILE_NAME,'Sigma_Theta',ierr)
          IF(ierr==1)THEN
             Sigma_Theta = 10.0_SP

      if (myid.eq.0) WRITE(3,'(A25,A50)')'Sigma_Theta', 'NOT DEFINED, USE 10.0'

          ENDIF ! end ierr

         ENDIF ! 1D or not



      if (myid.eq.0) WRITE(3,'(A12,I12)')'Nfreq=  ', Nfreq
      if (myid.eq.0) WRITE(3,'(A12,I12)')'Ntheta=  ', Ntheta
      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'FreqPeak=  ', FreqPeak
      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'FreqMin =  ', FreqMin
      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'FreqMax =  ', FreqMax
      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'Hmo     =  ', Hmo
      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'GammaTMA=  ', GammaTMA
      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'ThetaPeak= ', ThetaPeak
      if (myid.eq.0) WRITE(3,'(A13,F12.2)')'Sigma_Theta=', Sigma_Theta



      ENDIF ! end if data or tma
     ENDIF ! end absorbing-generating or left bc wave maker


      if (myid.eq.0) WRITE(3,*)'                                         '
      if (myid.eq.0) WRITE(3,*)'-------------- PERIODIC BC -----------------'



! south-north periodic boundary condition
      CALL READ_LOGICAL(PERIODIC,FILE_NAME,'PERIODIC',ierr)
      IF(ierr==1)THEN
        PERIODIC = .FALSE.
      ENDIF


      if (myid.eq.0) WRITE(3,'(A11,L2)')'PERIODIC:', PERIODIC



      if (myid.eq.0) WRITE(3,*)'                                         '
      if (myid.eq.0) WRITE(3,*)'-------------- SPONGE -----------------'


      CALL READ_LOGICAL(DIFFUSION_SPONGE,FILE_NAME,'DIFFUSION_SPONGE',ierr)
      IF(ierr==1)THEN
        DIFFUSION_SPONGE = .FALSE.
      ENDIF

      CALL READ_LOGICAL(DIRECT_SPONGE,FILE_NAME,'DIRECT_SPONGE',ierr)
      IF(ierr==1)THEN
        DIRECT_SPONGE = .FALSE.
      ENDIF

      CALL READ_LOGICAL(FRICTION_SPONGE,FILE_NAME,'FRICTION_SPONGE',ierr)
      IF(ierr==1)THEN
        FRICTION_SPONGE = .FALSE.
      ENDIF

      IF(DIRECT_SPONGE)THEN


      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'DIRECT_SPONGE IS USED'
         WRITE(3,'(A40)')'DIRECT_SPONGE IS USED'
      endif

       ENDIF

      IF(DIFFUSION_SPONGE)THEN


      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'DIFFUSION_SPONGE IS USED'
         WRITE(3,'(A40)')'DIFFUSION_SPONGE IS USED'
      endif


        CALL READ_FLOAT(Csp,FILE_NAME,'Csp',ierr)
      IF(ierr==1)THEN
        Csp = 0.1_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Csp Default:  0.1'
         WRITE(3,'(A40)')'Csp Default:  0.1'
      endif

       ENDIF


        if (myid.eq.0) WRITE(3,'(A22,F12.2)')'DIFFUSION_SPONGE Csp=', Csp

      ENDIF ! end diffusion_sponge

      IF(FRICTION_SPONGE)THEN


      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'FRICTION_SPONGE IS USED'
         WRITE(3,'(A40)')'FRICTION_SPONGE IS USED'
      endif


        CALL READ_FLOAT(CDsponge,FILE_NAME,'CDsponge',ierr)
      IF(ierr==1)THEN
        CDsponge = 5.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'CDsponge Default:  5.0'
         WRITE(3,'(A40)')'CDsponge Default:  5.0'
      endif

       ENDIF


        if (myid.eq.0) WRITE(3,'(A26,F12.2)')'FRICTION_SPONGE CDsponge=', CDsponge

      ENDIF  ! endif friction_sponge

      IF(DIFFUSION_SPONGE.OR.DIRECT_SPONGE.OR.FRICTION_SPONGE)THEN
        CALL READ_FLOAT(Sponge_west_width,FILE_NAME,'Sponge_west_width',ierr)
      IF(ierr==1)THEN
        Sponge_west_width = 0.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Sponge_west_width Default:  0.0'
         WRITE(3,'(A40)')'Sponge_west_width Default:  0.0'
      endif

       ENDIF

        CALL READ_FLOAT(Sponge_east_width,FILE_NAME,'Sponge_east_width',ierr)
      IF(ierr==1)THEN
        Sponge_east_width = 0.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Sponge_east_width Default:  0.0'
         WRITE(3,'(A40)')'Sponge_east_width Default:  0.0'
      endif

       ENDIF

        CALL READ_FLOAT(Sponge_south_width,FILE_NAME,'Sponge_south_width',ierr)
      IF(ierr==1)THEN
        Sponge_south_width = 0.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Sponge_south_width Default:  0.0'
         WRITE(3,'(A40)')'Sponge_south_width Default:  0.0'
      endif

       ENDIF

        CALL READ_FLOAT(Sponge_north_width,FILE_NAME,'Sponge_north_width',ierr)
      IF(ierr==1)THEN
        Sponge_north_width = 0.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Sponge_north_width Default:  0.0'
         WRITE(3,'(A40)')'Sponge_north_width Default:  0.0'
      endif

       ENDIF

        CALL READ_FLOAT(R_sponge,FILE_NAME,'R_sponge',ierr)
      IF(ierr==1)THEN
        R_sponge = 0.85_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'R_sponge Default:  0.85'
         WRITE(3,'(A40)')'R_sponge Default:  0.85'
      endif

       ENDIF

        CALL READ_FLOAT(A_sponge,FILE_NAME,'A_sponge',ierr)
      IF(ierr==1)THEN
        A_sponge = 5.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'A_sponge Default:  5.0'
         WRITE(3,'(A40)')'A_sponge Default:  5.0'
      endif

       ENDIF


        if (myid.eq.0) WRITE(3,'(A20,F12.2)')'Sponge_west_width =', Sponge_west_width
        if (myid.eq.0) WRITE(3,'(A20,F12.2)')'Sponge_east_width =', Sponge_east_width
        if (myid.eq.0) WRITE(3,'(A20,F12.2)')'Sponge_south_width=', Sponge_south_width
        if (myid.eq.0) WRITE(3,'(A20,F12.2)')'Sponge_north_width=', Sponge_north_width
        if (myid.eq.0) WRITE(3,'(A20,F12.2)')'R_sponge          =', R_sponge
        if (myid.eq.0) WRITE(3,'(A20,F12.2)')'A_sponge          =', A_sponge

       ENDIF ! endif sponge

! to avoid longshore current caused by extra momentum flux
! we can add bottom friction make momentum balance

      CALL READ_FLOAT(WaveMakerCd,FILE_NAME,'WaveMakerCd',ierr)
      IF(ierr==1)THEN
        WaveMakerCurrentBalance=.FALSE.

      if (myid.eq.0) WRITE(3,'(A40)')'No WavemakerCurrentBalance'

      ELSE
        WaveMakerCurrentBalance=.TRUE.

      if (myid.eq.0) WRITE(3,'(A15,F6.2)')'WaveMakerCd:', WaveMakerCd

      ENDIF



      if (myid.eq.0) WRITE(3,*)'                                         '
      if (myid.eq.0) WRITE(3,*)'-------- OBSTACLE and BREAKWATER -----------------'


! obstacle structures
      CALL READ_STRING(OBSTACLE_FILE,FILE_NAME,'OBSTACLE_FILE',ierr)
      IF(ierr==1)THEN
        OBSTACLE=.FALSE.

      if (myid.eq.0) WRITE(3,'(A15,A5)')'OBSTACLE_FILE:', 'NO'

      ELSE
        OBSTACLE=.TRUE.

      if (myid.eq.0) WRITE(3,'(A15,A50)')'OBSTACLE_FILE:', OBSTACLE_FILE

      ENDIF

! breakwater
      CALL READ_STRING(BREAKWATER_FILE,FILE_NAME,'BREAKWATER_FILE',ierr)
      IF(ierr==1)THEN
        BREAKWATER=.FALSE.

      if (myid.eq.0) WRITE(3,'(A20,A5)')'BREAKWATER_FILE:', 'NO'

      ELSE
        BREAKWATER=.TRUE.

      if (myid.eq.0) WRITE(3,'(A20,A50)')'BREAKWATER_FILE:', BREAKWATER_FILE

      ENDIF

! breakwater reflection stength

          CALL READ_FLOAT(BreakWaterAbsorbCoef,FILE_NAME,'BreakWaterAbsorbCoef',ierr)
      IF(ierr==1)THEN
        BreakWaterAbsorbCoef = 10.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'BreakWaterAbsorbCoef Default:  10.0'
         WRITE(3,'(A40)')'BreakWaterAbsorbCoef Default:  10.0'
      endif

       ELSE

      if (myid.eq.0) THEN
         WRITE(*,'(A40,F6.2)')'BreakWaterAbsorbCoef:', BreakWaterAbsorbCoef
         WRITE(3,'(A40,F6.2)')'BreakWaterAbsorbCoef:', BreakWaterAbsorbCoef
      endif

       ENDIF



      if (myid.eq.0) WRITE(3,*)'                                         '
      if (myid.eq.0) WRITE(3,*)'-------------- PHYSICS -----------------'


! physics
          CALL READ_LOGICAL(DISPERSION,FILE_NAME,'DISPERSION',ierr)
      IF(ierr==1)THEN
        DISPERSION = .TRUE.

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'DISPERSION Default:  DISPERSION'
         WRITE(3,'(A40)')'DISPERSION Default:  DISPERSION'
      endif

       ENDIF

          CALL READ_FLOAT(Gamma1,FILE_NAME,'Gamma1',ierr)
      IF(ierr==1)THEN
        Gamma1 = 1.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Gamma1 Default:  1.0: DISPERSION'
         WRITE(3,'(A40)')'Gamma1 Default:  1.0: DISPERSION'
      endif

       ENDIF


          CALL READ_FLOAT(Gamma2,FILE_NAME,'Gamma2',ierr)
      IF(ierr==1)THEN
        Gamma2 = 1.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'Gamma2 Default:  1.0: Full nonlinear'
         WRITE(3,'(A50)')'Gamma2 Default:  1.0: Full nonlinear'
      endif

       ENDIF

          CALL READ_FLOAT(Beta_ref,FILE_NAME,'Beta_ref',ierr)
      IF(ierr==1)THEN
        Beta_ref = - 0.531_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Beta_ref Default:  -0.531'
         WRITE(3,'(A40)')'Beta_ref Default:  -0.531'
      endif

       ENDIF



          CALL READ_FLOAT(Gamma3,FILE_NAME,'Gamma3',ierr)
      IF(ierr==1)THEN
        Gamma3 = 1.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'Gamma3 Default:  1.0: NOT fully linear'
         WRITE(3,'(A60)')'Gamma3 Default:  1.0: NOT fully linear'
      endif

       ENDIF


      if (myid.eq.0) WRITE(3,'(A20)')'Summary of Physics'




       if (myid.eq.0) WRITE(3,'(A10,F12.2)')'Gamma1 = ', Gamma1

       if (myid.eq.0) WRITE(3,'(A10,F12.2)')'Gamma2 = ', Gamma2
       if (myid.eq.0) WRITE(3,'(A10,F12.2)')'Beta_ref= ', Beta_ref
       if (myid.eq.0) WRITE(3,'(A10,F12.2)')'Gamma3 = ', Gamma3


      CALL READ_LOGICAL(VISCOSITY_BREAKING,FILE_NAME,'VISCOSITY_BREAKING',ierr)
      IF(ierr==1)THEN
        VISCOSITY_BREAKING = .TRUE.

      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'VISCOSITY_BREAKING Default:  VIS Breaking'
         WRITE(3,'(A60)')'VISCOSITY_BREAKING Default:  VIS Breaking'
      endif

       ENDIF

      IF(ROLLER) VISCOSITY_BREAKING = .TRUE.

      IF(VISCOSITY_BREAKING)THEN

       if (myid.eq.0) WRITE(3,*)'VISCOSITY_BREAKING IS USED'

      ENDIF

      CALL READ_FLOAT(SWE_ETA_DEP,FILE_NAME,'SWE_ETA_DEP',ierr)
      IF(ierr==1)THEN
        SWE_ETA_DEP = 0.80_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'SWE_ETA_DEP Default:  0.8'
         WRITE(3,'(A40)')'SWE_ETA_DEP Default:  0.8'
      endif

       ENDIF


      IF(VISCOSITY_BREAKING)THEN
        ! say nothing
      ELSE

       if (myid.eq.0) WRITE(3,'(A13,F12.2)')'SWE_ETA_DEP=', SWE_ETA_DEP

      ENDIF

      CALL READ_LOGICAL(IN_Cd,FILE_NAME,'FRICTION_MATRIX',ierr)
      IF(ierr==1)THEN
        IN_Cd = .FALSE.

      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'Friction_Matrix Default:  constant Cd'
         WRITE(3,'(A50)')'Friction_Matrix Default:  constant Cd'
      endif

       ENDIF

    IF(IN_Cd)THEN
        CALL READ_STRING(CD_FILE,FILE_NAME,'FRICTION_FILE',ierr) 
      IF(ierr==1)THEN

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'FRICTION_FILE NOT FOUND, Please specify Cd_file in input.txt'
         WRITE(3,'(A80)')'FRICTION_FILE NOT FOUND, Please specify Cd_file in input.txt'
      endif
       call MPI_FINALIZE ( ier )

        STOP
      ELSE

       if (myid.eq.0) WRITE(3,'(A15,A50)')'CD_FILE:', CD_FILE

      ENDIF

    ENDIF ! endif IN_Cd

 
      CALL READ_FLOAT(Cd_fixed,FILE_NAME,'Cd',ierr)
      IF(ierr==1)THEN
        Cd_fixed = 0.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'Cd_fixed Default:  0.0, possibly you used FRICTION_MATRIX'
         WRITE(3,'(A80)')'Cd_fixed Default:  0.0, possibly you used FRICTION_MATRIX'
      endif

       ENDIF
     

       if (myid.eq.0) WRITE(3,'(A35,F12.2)')'Cd_fixed (if you used fixed Cd) =', Cd_fixed



      if (myid.eq.0) WRITE(3,*)'                                         '
      if (myid.eq.0) WRITE(3,*)'-------------- NUMERICS -----------------'


! numerics schemes
      CALL READ_STRING(Time_Scheme,FILE_NAME,'Time_Scheme',ierr)
      IF(ierr==1)THEN
        Time_Scheme = 'Runge_Kutta'

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Time_Scheme Default:  Runge_Kutta'
         WRITE(3,'(A40)')'Time_Scheme Default:  Runge_Kutta'
      endif

       ENDIF


      if (myid.eq.0) WRITE(3,'(A13,A50)')'TIME_SCHEME:', TIME_SCHEME


      CALL READ_STRING(CONSTR,FILE_NAME,'CONSTRUCTION',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0) WRITE(3,'(A14,A50)')'CONSTRUCTION', 'NOT DEFINED, USE HLL'

        CONSTR='HLLC'
      ENDIF


      if (myid.eq.0) WRITE(3,'(A14,A50)')'CONSTRUCTION:', CONSTR


      CALL READ_STRING(HIGH_ORDER,FILE_NAME,'HIGH_ORDER',ierr)
      IF(ierr==1)THEN

      if (myid.eq.0)then
        WRITE(*,'(A12,A50)')'HIGH_ORDER', 'NOT DEFINED, USE FOURTH-ORDER'
        WRITE(3,'(A12,A50)')'HIGH_ORDER', 'NOT DEFINED, USE FOURTH-ORDER'
      endif

        HIGH_ORDER='FOURTH'        
      ENDIF


      if (myid.eq.0) WRITE(3,'(A12,A50)')'HIGH_ORDER:', HIGH_ORDER

! CFL
      CALL READ_FLOAT(CFL,FILE_NAME,'CFL',ierr)
      IF(ierr==1)THEN
        CFL = 0.5_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'CFL Default:  0.5'
         WRITE(3,'(A40)')'CFL Default:  0.5'
      endif

       ENDIF


      if (myid.eq.0) WRITE(3,'(A5,F12.2)')'CFL=', CFL


! DT_fixed
      CALL READ_FLOAT(DT_fixed,FILE_NAME,'DT_fixed',ierr)
      IF (ierr.ne.1) THEN
      FIXED_DT = .TRUE.

      if (myid.eq.0) then
           WRITE(3,'(A80)') 'use fixed DT, but judged by CFL. IF not satisfy CLF, DT/2...'
           WRITE(3,'(A12,F12.2)')'DT_fixed= ', DT_fixed
      endif

      ENDIF ! ierr\=1
  
! Froude Number Cap
      CALL READ_FLOAT(FroudeCap,FILE_NAME,'FroudeCap',ierr)
      IF(ierr==1)THEN
        FroudeCap = 3.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'FroudeCap Default:  3.0'
         WRITE(3,'(A40)')'FroudeCap Default:  3.0'
      endif

       ENDIF


      if (myid.eq.0) WRITE(3,'(A12,F12.2)')'FroudeCap=', FroudeCap


! MinDepth etc
      CALL READ_FLOAT(MinDepth,FILE_NAME,'MinDepth',ierr)
      IF(ierr==1)THEN
        MinDepth = 0.1_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'MinDepth Default:  0.1 m'
         WRITE(3,'(A40)')'MinDepth Default:  0.1 m'
      endif

       ENDIF

      CALL READ_FLOAT(MinDepthFrc,FILE_NAME,'MinDepthFrc',ierr)
      IF(ierr==1)THEN
        MinDepthFrc = 0.1_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'MinDepthFrc Default:  0.1 m'
         WRITE(3,'(A40)')'MinDepthFrc Default:  0.1 m'
      endif

       ENDIF

!  merge two parameters into the minimum one, change to min according to Harris 11/13/2023
       MinDepthFrc=MIN(MinDepthFrc,MinDepth)
       MinDepth=MinDepthFrc


      if (myid.eq.0) WRITE(3,'(A40)')'USE MIN(MinDepthFrc, MinDepth)'




      if (myid.eq.0) WRITE(3,'(A10,F12.6)')'MinDepth=', MinDepth



      if (myid.eq.0) WRITE(3,'(A13,F12.6)')'MinDepthFrc=', MinDepthFrc


! Lauren - Arrival Time - Wave height threshold (in m) to pick up arrival time

      CALL READ_LOGICAL(OUT_Time,FILE_NAME,'OUT_Time',ierr)
      IF(ierr==1)THEN
        OUT_Time = .FALSE.

      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'Dont record wave arrival time'
         WRITE(3,'(A60)')'Dont record wave arrival time'
      endif

       ELSE

      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'Record wave arrival time'
         WRITE(3,'(A60)')'Record wave arrival time'
      endif

       ENDIF

     IF(OUT_Time)THEN
      CALL READ_FLOAT(ArrTimeMin,FILE_NAME,'ArrTimeMinH',ierr)
      IF(ierr==1)THEN
        ArrTimeMin = 0.001 ! set equal to 0.1 cm if no threshold is given

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'ArrTimeMinH Default:  0.001 m'
         WRITE(3,'(A40)')'ArrTimeMinH Default:  0.001 m'
      endif

       ENDIF
       

      if (myid.eq.0) WRITE(3,'(A10,F12.6)')'ArrTimeMinH', ArrTimeMin

     ENDIF
 
! end Laurens modification


      if (myid.eq.0) WRITE(3,*)'                                         '
      if (myid.eq.0) WRITE(3,*)'--------- WAVE BREAKING -----------------'


! roller

      CALL READ_LOGICAL(ROLLER,FILE_NAME,'ROLLER_EFFECT',ierr)
      IF(ierr==1)THEN
        ROLLER = .FALSE.
      ENDIF

       IF(ROLLER)THEN
         ROLLER_SWITCH = 1.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'ROLLER_EFFECT:  INCLUDED'
         WRITE(3,'(A40)')'ROLLER_EFFECT:  INCLUDED'
      endif

       ELSE
         ROLLER_SWITCH = ZERO

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'ROLLER_EFFECT:  NO'
         WRITE(3,'(A40)')'ROLLER_EFFECT:  NO'
      endif

       ENDIF

! end roller

! show breaking
      CALL READ_LOGICAL(SHOW_BREAKING,FILE_NAME,'SHOW_BREAKING',ierr)
      IF(ierr==1)THEN
        SHOW_BREAKING = .TRUE.

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'SHOW_BREAKING Default:  TRUE'
         WRITE(3,'(A40)')'SHOW_BREAKING Default:  TRUE'
      endif

       ENDIF
	
      IF(VISCOSITY_BREAKING) SHOW_BREAKING = .TRUE.

      IF(SHOW_BREAKING)THEN
      CALL READ_FLOAT(Cbrk1,FILE_NAME,'Cbrk1',ierr)
      IF(ierr==1)THEN
        Cbrk1 = 0.65_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Cbrk1 Default:  0.65'
         WRITE(3,'(A40)')'Cbrk1 Default:  0.65'
      endif

       ENDIF

      IF(VISCOSITY_BREAKING)THEN

      if (myid.eq.0) WRITE(3,'(A8,F12.6)')'Cbrk1 =', Cbrk1

      ENDIF

      CALL READ_FLOAT(Cbrk2,FILE_NAME,'Cbrk2',ierr)
      IF(ierr==1)THEN
        Cbrk2 = 0.35_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Cbrk2 Default:  0.35'
         WRITE(3,'(A40)')'Cbrk2 Default:  0.35'
      endif

       ENDIF

      IF(VISCOSITY_BREAKING)THEN

      if (myid.eq.0) WRITE(3,'(A8,F12.6)')'Cbrk2 =', Cbrk2

      ENDIF

      CALL READ_FLOAT(WAVEMAKER_Cbrk,FILE_NAME,'WAVEMAKER_Cbrk',ierr)
      IF(ierr==1)THEN
        WAVEMAKER_Cbrk = 1.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'WAVEMAKER_Cbrk Default:  1.0'
         WRITE(3,'(A40)')'WAVEMAKER_Cbrk Default:  1.0'
      endif

       ENDIF


      if (myid.eq.0) WRITE(3,'(A18,F17.6)')'WAVEMAKER_Cbrk =', WAVEMAKER_Cbrk

      ENDIF

	![ykchoi(08.18.2015) : for viscosity of wavemaker
      CALL READ_LOGICAL(WAVEMAKER_VIS,FILE_NAME,'WAVEMAKER_VIS',ierr)  
      IF(ierr==1)THEN
        WAVEMAKER_VIS = .FALSE.

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'WAVEMAKER_VIS Default:  FALSE'
         WRITE(3,'(A40)')'WAVEMAKER_VIS Default:  FALSE'
      endif

       ENDIF
	
	IF( VISCOSITY_BREAKING .AND. WAVEMAKER_VIS ) THEN

	  IF (myid.eq.0) then
	     WRITE(*,*) "==============================================="
	     WRITE(*,*)  "STOP :: VISCOSITY_BREAKING=T, WAVEMAKER_VIS=T"
	     WRITE(*,*) "==============================================="
          ENDIF
          call MPI_FINALIZE ( ier )

      ENDIF

      IF(WAVEMAKER_VIS)THEN

       if (myid.eq.0) WRITE(3,*)'WAVEMAKER_VIS'

      ENDIF

!  suggest dont use wavemaker_vis 04/30
      IF(WAVEMAKER_VIS)THEN
      CALL READ_FLOAT(visbrk,FILE_NAME,'visbrk',ierr)
	CALL READ_FLOAT(WAVEMAKER_visbrk,FILE_NAME,'WAVEMAKER_visbrk',ierr)

      if (myid.eq.0) WRITE(3,'(A14,F12.6)')'visbrk =', visbrk
      if (myid.eq.0) WRITE(3,'(A8,F12.6)')'WAVEMAKER_visbrk =', WAVEMAKER_visbrk

	ENDIF


      if (myid.eq.0) WRITE(3,*)'                                         '
      if (myid.eq.0) WRITE(3,*)'-------- WAVE-AVERAGED PROPERTY -----------------'


      CALL READ_FLOAT(T_INTV_mean,FILE_NAME,'T_INTV_mean',ierr)
      IF(ierr==1)THEN
        T_INTV_mean = LARGE

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'T_INTV_mean Default:  LARGE'
         WRITE(3,'(A40)')'T_INTV_mean Default:  LARGE'
      endif

       ENDIF

	CALL READ_FLOAT(STEADY_TIME,FILE_NAME,'STEADY_TIME',ierr)
      IF(ierr==1)THEN
        STEADY_TIME = LARGE

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'STEADY_TIME Default:  LARGE'
         WRITE(3,'(A40)')'STEADY_TIME Default:  LARGE'
      endif

       ENDIF

      CALL READ_FLOAT(C_smg,FILE_NAME,'C_smg',ierr)
      IF(ierr==1)THEN
        C_smg = 0.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'C_smg Default:  0.0'
         WRITE(3,'(A40)')'C_smg Default:  0.0'
      endif

       ENDIF


      if (myid.eq.0) WRITE(3,'(A14,F12.6)')'T_INTV_mean =', T_INTV_mean
      if (myid.eq.0) WRITE(3,'(A14,F12.6)')'STEADY_TIME =', STEADY_TIME
      if (myid.eq.0) WRITE(3,'(A8,F12.6)')'C_smg =', C_smg

	
      CALL READ_FLOAT(nu_bkg,FILE_NAME,'nu_bkg',ierr)
      IF(ierr==1)THEN
        nu_bkg = 0.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'nu_bkg Default:  0.0'
         WRITE(3,'(A40)')'nu_bkg Default:  0.0'
      endif

       ENDIF


  ! end coupling file



      if (myid.eq.0) WRITE(3,*)'                                         '
      if (myid.eq.0) WRITE(3,*)'-------------- OUTPUT INFO -----------------'


! result folder
      CALL READ_STRING(RESULT_FOLDER,FILE_NAME,'RESULT_FOLDER',ierr)
      IF(ierr==1)THEN
        RESULT_FOLDER = './output/'
      ENDIF


      if (myid.eq.0) WRITE(3,'(A15,A50)')'RESULT_FOLDER:', RESULT_FOLDER


![mayhl 17/06/12
! Note 1: Serial code only output is ASCII
! Note 2: Placing code after result folder creation seems to delay
!         writing to LOG.txt till after simulation is completed.

      CALL READ_STRING(FIELD_IO_TYPE,FILE_NAME,'FIELD_IO_TYPE',ierr)

      IF(ierr.EQ.1) FIELD_IO_TYPE = 'ASCII'
      
      IF (myid.EQ.0) WRITE(3,*) 'FIELD_IO_TYPE = ' , FIELD_IO_TYPE

! mayhl]
 


! create result folder
      MKFOLDER = "mkdir -p "//TRIM(RESULT_FOLDER)

      IF (myid.eq.0) THEN

        CALL SYSTEM(TRIM(MKFOLDER))

      ENDIF


! station files
      CALL READ_INTEGER(NumberStations,FILE_NAME,'NumberStations',ierr)
      IF(NumberStations>0)THEN
      CALL READ_STRING(STATIONS_FILE,FILE_NAME,'STATIONS_FILE',ierr)
      ENDIF

! output parameters
      CALL READ_INTEGER(OUTPUT_RES,FILE_NAME,'OUTPUT_RES',ierr)
      IF(ierr==1)THEN
        OUTPUT_RES = 1

      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'OUTPUT_RES NOT FOUND, OUTPUT_RES=1: full resolution'
         WRITE(3,'(A60)')'OUTPUT_RES NOT FOUND, OUTPUT_RES=1: full resolution'
      endif

       ENDIF


       if (myid.eq.0) WRITE(3,'(A15,I10)')'OUTPUT_RES',OUTPUT_RES


      CALL READ_LOGICAL(OUT_DEPTH,FILE_NAME,'DEPTH_OUT',ierr)
      CALL READ_LOGICAL(OUT_U,FILE_NAME,'U',ierr)
      CALL READ_LOGICAL(OUT_V,FILE_NAME,'V',ierr)
      CALL READ_LOGICAL(OUT_ETA,FILE_NAME,'ETA',ierr)

      CALL READ_LOGICAL(OUT_Hmax,FILE_NAME,'Hmax',ierr)
      CALL READ_LOGICAL(OUT_Hmin,FILE_NAME,'Hmin',ierr)
      CALL READ_LOGICAL(OUT_Umax,FILE_NAME,'Umax',ierr)
      CALL READ_LOGICAL(OUT_MFmax,FILE_NAME,'MFmax',ierr)
      CALL READ_LOGICAL(OUT_VORmax,FILE_NAME,'VORmax',ierr)
      CALL READ_LOGICAL(OUT_MASK,FILE_NAME,'MASK',ierr)
      CALL READ_LOGICAL(OUT_MASK9,FILE_NAME,'MASK9',ierr)
      CALL READ_LOGICAL(OUT_Umean,FILE_NAME,'Umean',ierr)
      CALL READ_LOGICAL(OUT_Vmean,FILE_NAME,'Vmean',ierr)
      CALL READ_LOGICAL(OUT_ETAmean,FILE_NAME,'ETAmean',ierr)
      CALL READ_LOGICAL(OUT_WaveHeight,FILE_NAME,'WaveHeight',ierr)
      CALL READ_LOGICAL(OUT_SXL,FILE_NAME,'SXL',ierr)
      CALL READ_LOGICAL(OUT_SXR,FILE_NAME,'SXR',ierr)
      CALL READ_LOGICAL(OUT_SYL,FILE_NAME,'SYL',ierr)
      CALL READ_LOGICAL(OUT_SYR,FILE_NAME,'SYR',ierr)
      CALL READ_LOGICAL(OUT_SourceX,FILE_NAME,'SourceX',ierr)
      CALL READ_LOGICAL(OUT_SourceY,FILE_NAME,'SourceY',ierr)
      CALL READ_LOGICAL(OUT_FrcX,FILE_NAME,'FrcX',ierr)
      CALL READ_LOGICAL(OUT_FrcY,FILE_NAME,'FrcY',ierr)
      CALL READ_LOGICAL(OUT_BrkdisX,FILE_NAME,'BrkdisX',ierr)
      CALL READ_LOGICAL(OUT_BrkdisY,FILE_NAME,'BrkdisY',ierr)
      CALL READ_LOGICAL(OUT_P,FILE_NAME,'P',ierr)
      CALL READ_LOGICAL(OUT_Q,FILE_NAME,'Q',ierr)
      CALL READ_LOGICAL(OUT_Fx,FILE_NAME,'Fx',ierr)
      CALL READ_LOGICAL(OUT_Fy,FILE_NAME,'Fy',ierr)
      CALL READ_LOGICAL(OUT_Gx,FILE_NAME,'Gx',ierr)
      CALL READ_LOGICAL(OUT_Gy,FILE_NAME,'Gy',ierr)
      CALL READ_LOGICAL(OUT_AGE,FILE_NAME,'AGE',ierr)
      CALL READ_LOGICAL(OUT_ROLLER,FILE_NAME,'ROLLER',ierr)
      CALL READ_LOGICAL(OUT_UNDERTOW,FILE_NAME,'UNDERTOW',ierr)
      CALL READ_LOGICAL(OUT_NU,FILE_NAME,'OUT_NU',ierr)
      CALL READ_LOGICAL(OUT_TMP,FILE_NAME,'TMP',ierr) 
      CALL READ_LOGICAL(OUT_Radiation,FILE_NAME,'Radiation',ierr) 

!ykchoi
!	CALL READ_FLOAT(EtaBlowVal,FILE_NAME,'EtaBlowVal',ierr)
      IF(ierr==1)THEN
        EtaBlowVal = 10.0_SP

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'EtaBlowVal Default:  100xmax_depth'
         WRITE(3,'(A40)')'EtaBlowVal Default:  100xmax_depth'
      endif

       ENDIF

!  fyshi set blowup value is 100xmax_depth in init.F  


      if (myid.eq.0)   then

      WRITE(3,'(A15,L2)')'OUT_DEPTH',OUT_DEPTH
      WRITE(3,'(A15,L2)')'OUT_U',OUT_U
      WRITE(3,'(A15,L2)')'OUT_V',OUT_V
      WRITE(3,'(A15,L2)')'OUT_ETA',OUT_ETA

      WRITE(3,'(A15,L2)')'OUT_Hmax',OUT_Hmax
      WRITE(3,'(A15,L2)')'OUT_Hmin',OUT_Hmin
      WRITE(3,'(A15,L2)')'OUT_Umax',OUT_Umax
      WRITE(3,'(A15,L2)')'OUT_MFmax',OUT_MFmax
      WRITE(3,'(A15,L2)')'OUT_VORmax',OUT_VORmax
      WRITE(3,'(A15,L2)')'OUT_MASK',OUT_MASK
      WRITE(3,'(A15,L2)')'OUT_MASK9',OUT_MASK9
      WRITE(3,'(A15,L2)')'OUT_Umean',OUT_Umean
      WRITE(3,'(A15,L2)')'OUT_Vmean',OUT_Vmean
      WRITE(3,'(A15,L2)')'OUT_ETAmean',OUT_ETAmean
      WRITE(3,'(A15,L2)')'OUT_WaveHeight',OUT_WaveHeight
      WRITE(3,'(A15,L2)')'OUT_SXL',OUT_SXL
      WRITE(3,'(A15,L2)')'OUT_SXR',OUT_SXR
      WRITE(3,'(A15,L2)')'OUT_SYL',OUT_SYL
      WRITE(3,'(A15,L2)')'OUT_SYR',OUT_SYR
      WRITE(3,'(A15,L2)')'OUT_SourceX',OUT_SourceX
      WRITE(3,'(A15,L2)')'OUT_SourceY',OUT_SourceY
      WRITE(3,'(A15,L2)')'OUT_P',OUT_P
      WRITE(3,'(A15,L2)')'OUT_Q',OUT_Q
      WRITE(3,'(A15,L2)')'OUT_Fx',OUT_Fx
      WRITE(3,'(A15,L2)')'OUT_Fy',OUT_Fy
      WRITE(3,'(A15,L2)')'OUT_Gx',OUT_Gx
      WRITE(3,'(A15,L2)')'OUT_Gy',OUT_Gy
      WRITE(3,'(A15,L2)')'OUT_AGE',OUT_AGE
      WRITE(3,'(A15,L2)')'OUT_ROLLER',OUT_ROLLER
      WRITE(3,'(A15,L2)')'OUT_UNDERTOW',OUT_UNDERTOW
      WRITE(3,'(A15,L2)')'OUT_NU',OUT_NU
      WRITE(3,'(A15,L2)')'OUT_TMP',OUT_TMP
      WRITE(3,'(A15,L2)')'OUT_TIME',OUT_Time

      endif




      if (myid.eq.0) WRITE(3,*)'                                         '
      if (myid.eq.0) WRITE(3,*)' --------------input end --------------' 
      if (myid.eq.0) WRITE(3,*)'                                         '


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!START!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Oct, 2021
!Salatin, R., Chen, Q., Bak, A. S., Shi, F., & Brandt, S. R. (2021). Effects of
!wave coherence on longshore variability of nearshore wave processes. Journal
!of Geophysical Research: Oceans,126, e2021JC017641.
!https://doi.org/10.1029/2021JC017641

    IF(WaveMaker(1:10)=='WK_NEW_IRR')THEN
        CALL READ_FLOAT(Xc_WK,FILE_NAME,'Xc_WK',ierr)
        IF(ierr==1)THEN

            if (myid.eq.0) THEN
                WRITE(*,'(A60)')'Xc_WK NOT FOUND, specify Xc_WK in input.txt'
                WRITE(3,'(A60)')'Xc_WK NOT FOUND, specify Xc_WK in input.txt'
            endif
            call MPI_FINALIZE ( ier )

            STOP
        ENDIF

        CALL READ_FLOAT(Yc_WK,FILE_NAME,'Yc_WK',ierr)

        IF(ierr==1)THEN
            Yc_WK = ZERO

            if (myid.eq.0) THEN
                WRITE(*,'(A40)')'Yc_WK defalt: 0.0'
                WRITE(3,'(A40)')'Yc_WK defalt: 0.0'
            endif

        ENDIF

        CALL READ_FLOAT(DEP_WK,FILE_NAME,'DEP_WK',ierr)

        IF(ierr==1)THEN

            if (myid.eq.0) THEN
                WRITE(*,'(A60)')'DEP_WK NOT FOUND, specify DEP_WK in input.txt'
                WRITE(3,'(A60)')'DEP_WK NOT FOUND, specify DEP_WK in input.txt'
            endif
            call MPI_FINALIZE ( ier )

            STOP
        ENDIF

        CALL READ_FLOAT(Time_ramp,FILE_NAME,'Time_ramp',ierr)

        IF(ierr==1)THEN
            Time_ramp = 0.0_SP

            if (myid.eq.0) THEN
                WRITE(*,'(A40)')'Time_ramp Default:  0.0'
                WRITE(3,'(A40)')'Time_ramp Default:  0.0'
            endif

        ENDIF

        CALL READ_FLOAT(Delta_WK,FILE_NAME,'Delta_WK',ierr)

        IF(ierr==1)THEN
            Delta_WK = 0.5_SP

            if (myid.eq.0) THEN
                WRITE(*,'(A40)')'Delta_WK Default:  0.5'
                WRITE(3,'(A40)')'Delta_WK Default:  0.5'
            endif

        ENDIF

        CALL READ_FLOAT(FreqPeak,FILE_NAME,'FreqPeak',ierr)

        IF(ierr==1)THEN

            if (myid.eq.0) THEN
                WRITE(*,'(A80)')'FreqPeak NOT FOUND, specify FreqPeak in input.txt'
                WRITE(3,'(A80)')'FreqPeak NOT FOUND, specify FreqPeak in input.txt'
            endif
            call MPI_FINALIZE ( ier )

            STOP
        ENDIF

        CALL READ_FLOAT(FreqMin,FILE_NAME,'FreqMin',ierr)

        IF(ierr==1)THEN

            if (myid.eq.0) THEN
                WRITE(*,'(A80)')'FreqMin NOT FOUND, specify FreqMin in input.txt'
                WRITE(3,'(A80)')'FreqMin NOT FOUND, specify FreqMin in input.txt'
            endif
            call MPI_FINALIZE ( ier )

            STOP
        ENDIF

        CALL READ_FLOAT(FreqMax,FILE_NAME,'FreqMax',ierr)

        IF(ierr==1)THEN

            if (myid.eq.0) THEN
                WRITE(*,'(A80)')'FreqMax NOT FOUND, specify FreqMax in input.txt'
                WRITE(3,'(A80)')'FreqMax NOT FOUND, specify FreqMax in input.txt'
            endif
            call MPI_FINALIZE ( ier )

            STOP
        ENDIF

        CALL READ_FLOAT(Hmo,FILE_NAME,'Hmo',ierr)

        IF(ierr==1)THEN

            if (myid.eq.0) THEN
                WRITE(*,'(A50)')'Hmo NOT FOUND, specify Hmo in input.txt'
                WRITE(3,'(A50)')'Hmo NOT FOUND, specify Hmo in input.txt'
            endif
            call MPI_FINALIZE ( ier )

            STOP
        ENDIF

        CALL READ_FLOAT(GammaTMA,FILE_NAME,'GammaTMA',ierr)

        IF(ierr==1)THEN
            GammaTMA = 3.3

            if (myid.eq.0) WRITE(3,'(A14,A50)')'GammaTMA', 'NOT DEFINED, USE 3.3'

        ENDIF

        CALL READ_INTEGER(Nfreq,FILE_NAME,'Nfreq',ierr)

        IF(ierr==1)THEN
            Nfreq = 1080

            if (myid.eq.0) WRITE(3,'(A25,A50)')'Nfreq for TMA or JON', 'NOT DEFINED, USE 24*45=1080'

        ENDIF

        CALL READ_INTEGER(Ntheta,FILE_NAME,'Ntheta',ierr)

        IF(ierr==1)THEN
            Ntheta = 24

            if (myid.eq.0) WRITE(3,'(A25,A50)')'Ntheta for TMA or JON', 'NOT DEFINED, USE 24'

        ENDIF

        CALL READ_FLOAT(ThetaPeak,FILE_NAME,'ThetaPeak',ierr)

        IF(ierr==1)THEN
            ThetaPeak = 0.0_SP

            if (myid.eq.0) WRITE(3,'(A25,A50)')'ThetaPeak', 'NOT DEFINED, USE 0.0'

        ENDIF ! end ierr

        CALL READ_FLOAT(Sigma_Theta,FILE_NAME,'Sigma_Theta',ierr)

        IF(ierr==1)THEN
            Sigma_Theta = 10.0_SP

            if (myid.eq.0) WRITE(3,'(A25,A50)')'Sigma_Theta', 'NOT DEFINED, USE 10.0'

        ENDIF ! end ierr

        CALL READ_FLOAT(Ywidth_WK,FILE_NAME,'Ywidth_WK',ierr)

        IF(ierr==1)THEN
            Ywidth_WK = LARGE

            if (myid.eq.0) THEN
                WRITE(*,'(A50)')'Ywidth_WK Default:  LARGE'
                WRITE(3,'(A50)')'Ywidth_WK Default:  LARGE'
            endif

        ENDIF

        CALL READ_FLOAT(alpha_c,FILE_NAME,'alpha_c',ierr)

        IF(ierr==1)THEN
            alpha_c = 0.0_SP

            if (myid.eq.0) THEN
                WRITE(*,'(A50)')'Wave Coherence Percentage:  0.0%'
                WRITE(3,'(A50)')'Wave Coherence Percentage:  0.0%'
            endif

        ENDIF


        if (myid.eq.0) WRITE(3,'(A12,F12.2)')'Xc_WK   =  ', Xc_WK
        if (myid.eq.0) WRITE(3,'(A10,F12.2)')'Yc_WK   =', Yc_WK
        if (myid.eq.0) WRITE(3,'(A11,F12.2)')'Ywidth_WK=', Ywidth_WK
        if (myid.eq.0) WRITE(3,'(A12,F12.2)')'DEP_WK  =  ', DEP_WK
        if (myid.eq.0) WRITE(3,'(A12,F12.2)')'Time_ramp= ', Time_ramp
        if (myid.eq.0) WRITE(3,'(A12,F12.2)')'Delta_WK=  ', Delta_WK
        if (myid.eq.0) WRITE(3,'(A12,I12)')'Nfreq=  ', Nfreq
        if (myid.eq.0) WRITE(3,'(A12,I12)')'Ntheta=  ', Ntheta
        if (myid.eq.0) WRITE(3,'(A12,F12.2)')'FreqPeak=  ', FreqPeak
        if (myid.eq.0) WRITE(3,'(A12,F12.2)')'FreqMin =  ', FreqMin
        if (myid.eq.0) WRITE(3,'(A12,F12.2)')'FreqMax =  ', FreqMax
        if (myid.eq.0) WRITE(3,'(A12,F12.2)')'Hmo     =  ', Hmo
        if (myid.eq.0) WRITE(3,'(A12,F12.2)')'GammaTMA=  ', GammaTMA
        if (myid.eq.0) WRITE(3,'(A12,F12.2)')'ThetaPeak= ', ThetaPeak
        if (myid.eq.0) WRITE(3,'(A13,F12.2)')'Sigma_Theta=', Sigma_Theta
        if (myid.eq.0) WRITE(3,'(A13,F12.2)')'alpha_c=', alpha_c


    ENDIF ! IF(WaveMaker(1:10)=='WK_NEW_IRR')THEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE READ_INPUT


!-------------------------------------------------------------------------------------
!
!    STATIONS is a subroutine to write station data
! Fengyan Shi modified based on Jeff Harris for Spherical
! here simply specify grid number i and j instead of x and y
!
! HISTORY: 
!    09/16/2011  Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE STATIONS
     USE GLOBAL
     USE INPUT_READ

     IMPLICIT NONE

     INTEGER :: iunit
     INTEGER, DIMENSION(:),ALLOCATABLE :: BufferCount
     REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: ZUV_Buffer
     REAL(SP) :: dum1,dum2
     REAL(SP) :: eta_sta,u_sta,v_sta
     CHARACTER(LEN=80)::FILE_NAME=' '
     CHARACTER(LEN=80)::TMP_NAME=' '
     CHARACTER(LEN=80)::FDIR=' '
     LOGICAL :: FirstCallStation = .TRUE.
     SAVE FirstCallStation, BufferCount,ZUV_Buffer


! initialize stations
     FDIR=TRIM(RESULT_FOLDER)
     if (FirstCallStation) then
       FirstCallStation = .FALSE.
       ALLOCATE(ista(NumberStations),&
                jsta(NumberStations),&
                nsta(NumberStations))


       ALLOCATE(BufferCount(NumberStations), &
            ZUV_Buffer(StationOutputBuffer,NumberStations,4))


       BufferCount = 0
! calculate how many output components

! check existing

 INQUIRE(FILE=TRIM(STATIONS_FILE),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN

   IF(MYID==0)  &
   WRITE(*,*) TRIM(STATIONS_FILE), ' THE STATION FILE CANNOT BE FOUND. STOP'
   CALL MPI_FINALIZE (ier)
   STOP

  ENDIF  ! exist
              
       open(100,FILE=TRIM(STATIONS_FILE))
       do i=1,NumberStations
          read(100,*) dum1,dum2

![---ykchoi Jan/23/2018
          !ista(i) = Nghost+dum1-npx*Mglob/px
          !jsta(i) = Nghost+dum2-npy*Nglob/py
		ista(i) = Nghost+dum1-( iista - 1 )
		jsta(i) = Nghost+dum2-( jjsta - 1 )
!---ykchoi Jan/23/2018]
          if ((ista(i).ge.Ibeg).and.(ista(i).le.Iend).and.&
              (jsta(i).ge.Jbeg).and.(jsta(i).le.Jend)) then
             nsta(i) = 1
             write(file_name(1:4),'(I4.4)') i
             TMP_NAME = TRIM(FDIR)//'sta_'//TRIM(FILE_NAME)
             iunit=100+i
             open(iunit,FILE=TMP_NAME)
          else
             nsta(i) = 0
          endif


       enddo
     endif

! write to stations

     do i=1,NumberStations
       if (nsta(i).eq.1) then
          iunit=100+i
           IF(mask(ista(i),jsta(i))<1)THEN
             eta_sta=ZERO
             u_sta=ZERO
             v_sta=ZERO

           ELSE  ! to avoid topography on nested water surface  
             eta_sta=eta(ista(i),jsta(i))
             u_sta=u(ista(i),jsta(i))
             v_sta=v(ista(i),jsta(i)) 

           ENDIF

            IF(BufferCount(i)<StationOutputBuffer.AND.TIME<TOTAL_TIME)THEN
              BufferCount(i) = BufferCount(i) +1
              ZUV_BUFFER(BufferCount(i),i,1)=time
              ZUV_BUFFER(BufferCount(i),i,2)=eta_sta
              ZUV_BUFFER(BufferCount(i),i,3)=u_sta
              ZUV_BUFFER(BufferCount(i),i,4)=v_sta

            ELSE
              DO j=1,StationOutputBuffer

                write (iunit,'(E21.10E4, 3E16.5E4)') (ZUV_BUFFER(j,i,k),k=1,4)

              ENDDO              
              BufferCount(i) = 0
              ZUV_BUFFER(:,i,1)=TOTAL_TIME

              ZUV_BUFFER(:,i,2:4)=ZERO



      if (myid.eq.0) WRITE(*,*) 'WRITE OUT STATION ...'



            ENDIF ! end buffer

       endif
     enddo

! close station files
     if (TIME.ge.TOTAL_TIME) then
       do i=1,NumberStations
          if (nsta(i).eq.1) then
             iunit=100+i
             close(iunit)
          endif
       enddo
     endif

END SUBROUTINE STATIONS



!-------------------------------------------------------------------------------------
!
!    PREVIEW is subroutine for print-out of field data
!
!  HISTORY:
!    05/01/2010  Fengyan Shi
!    06/01/2015  Young-Kwang Choi, change file number to 5 digits, 
!                        such as eta_00001
!
!-------------------------------------------------------------------------------------
SUBROUTINE PREVIEW
     USE GLOBAL



     IMPLICIT NONE

     CHARACTER(LEN=80)::FILE_NAME=' '
     CHARACTER(LEN=80)::FILE_NAME_MEAN=' '
     CHARACTER(LEN=80)::TMP_NAME=' '
     CHARACTER(LEN=80)::FDIR=' '

     FDIR=TRIM(RESULT_FOLDER)

     ICOUNT=ICOUNT+1


        if (myid.eq.0)then
        WRITE(3,102)'PRINTING FILE NO.', icount, ' TIME/TOTAL: ', TIME,'/',Total_Time
        WRITE(*,102)'PRINTING FILE NO.', icount, ' TIME/TOTAL: ', TIME,'/',Total_Time        
        endif


102     FORMAT(A20,I6,A14,F12.3,A2,F12.3)

        !ykchoi
	  !itmp1=mod(icount/1000,10)
        !itmp2=mod(icount/100,10)
        !itmp3=mod(icount/10,10)
        !itmp4=mod(icount,10)
	   itmp1=mod(icount/10000,10)
	   itmp2=mod(icount/1000,10)
	   itmp3=mod(icount/100,10)
	   itmp4=mod(icount/10,10)
	   itmp5=mod(icount,10)

        write(file_name(1:1),'(I1)')itmp1
        write(file_name(2:2),'(I1)')itmp2
        write(file_name(3:3),'(I1)')itmp3
        write(file_name(4:4),'(I1)')itmp4
    	  write(file_name(5:5),'(I1)')itmp5   !ykchoi

     IF(ICOUNT==1)THEN
     IF(OUT_DEPTH.OR.BREAKWATER)THEN
        TMP_NAME = TRIM(FDIR)//'dep.out'
        call PutFile(TMP_NAME,DEPTH)
        TMP_NAME = TRIM(FDIR)//'cd_breakwater.out'
        call PutFile(TMP_NAME,CD_breakwater)
     ENDIF
     ENDIF
![ykchoi
     write(10000,*)time, dt
!ykchoi]

     IF(OUT_ETA)THEN
        TMP_NAME = TRIM(FDIR)//'eta_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,Eta)
     ENDIF




        TMP_NAME = TRIM(FDIR)//'Ax_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,grdAx)
        TMP_NAME = TRIM(FDIR)//'Ay_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,grdAy)
        TMP_NAME = TRIM(FDIR)//'Bx_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,grdBx)
        TMP_NAME = TRIM(FDIR)//'By_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,grdBy)


     IF(OUT_Hmax)THEN
        TMP_NAME = TRIM(FDIR)//'hmax_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,HeightMax)
     ENDIF

     IF(OUT_Hmin)THEN
        TMP_NAME = TRIM(FDIR)//'hmin_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,HeightMin)
     ENDIF

     IF(OUT_Umax)THEN
        TMP_NAME = TRIM(FDIR)//'umax_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,VelocityMax)
     ENDIF
     
     IF(OUT_MFmax)THEN                                                                                            
        TMP_NAME = TRIM(FDIR)//'MFmax_'//TRIM(FILE_NAME)                                                          
        call PutFile(TMP_NAME,MomentumFluxMax)                                                                              
     ENDIF      
     
     IF(OUT_VORmax)THEN                                                                                            
        TMP_NAME = TRIM(FDIR)//'VORmax_'//TRIM(FILE_NAME)                                                          
        call PutFile(TMP_NAME,VorticityMax)                                                                              
     ENDIF            
     
     IF(OUT_U)THEN
        TMP_NAME = TRIM(FDIR)//'u_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,U)
     ENDIF

     IF(OUT_V)THEN
        TMP_NAME = TRIM(FDIR)//'v_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,V)
     ENDIF

     IF(OUT_MASK)THEN
        TMP_NAME = TRIM(FDIR)//'mask_'//TRIM(FILE_NAME)
        Int2Flo=MASK
        call PutFile(TMP_NAME,Int2Flo)
     ENDIF

     IF(OUT_MASK9)THEN
        TMP_NAME = TRIM(FDIR)//'mask9_'//TRIM(FILE_NAME)
        Int2Flo=MASK9
        call PutFile(TMP_NAME,Int2Flo)
     ENDIF

210   FORMAT(5000I3)

     IF(OUT_P)THEN
        TMP_NAME = TRIM(FDIR)//'p_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,P(1:Mloc,1:Nloc))
     ENDIF

     IF(OUT_Q)THEN
        TMP_NAME = TRIM(FDIR)//'q_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,Q(1:Mloc,1:Nloc))
     ENDIF


     IF(OUT_AGE)THEN
      IF(SHOW_BREAKING)THEN
        TMP_NAME = TRIM(FDIR)//'age_'//TRIM(FILE_NAME)
         call PutFile(TMP_NAME,AGE_BREAKING)
      ENDIF
     ENDIF

     IF(OUT_ROLLER)THEN
        TMP_NAME = TRIM(FDIR)//'roller_'//TRIM(FILE_NAME)
         call PutFile(TMP_NAME,ROLLER_FLUX)
     ENDIF

     IF(OUT_UNDERTOW)THEN
        TMP_NAME = TRIM(FDIR)//'U_undertow_'//TRIM(FILE_NAME)
         call PutFile(TMP_NAME,UNDERTOW_U)
        TMP_NAME = TRIM(FDIR)//'V_undertow_'//TRIM(FILE_NAME)
         call PutFile(TMP_NAME,UNDERTOW_V)
     ENDIF

      IF(VISCOSITY_BREAKING)THEN
       IF(OUT_NU)THEN
        TMP_NAME = TRIM(FDIR)//'nubrk_'//TRIM(FILE_NAME)
         call PutFile(TMP_NAME,nu_break)
       ENDIF
!        TMP_NAME = TRIM(FDIR)//'etat_'//TRIM(FILE_NAME)
!         call PutFile(TMP_NAME,etat) 
      ENDIF

       IF(OUT_FrcX)THEN
        TMP_NAME = TRIM(FDIR)//'FrcInsX_'//TRIM(FILE_NAME)
         call PutFile(TMP_NAME,FrcInsX)
       ENDIF
       IF(OUT_FrcY)THEN
        TMP_NAME = TRIM(FDIR)//'FrcInsY_'//TRIM(FILE_NAME)
         call PutFile(TMP_NAME,FrcInsY)
       ENDIF
       IF(OUT_BrkdisX)THEN
        TMP_NAME = TRIM(FDIR)//'BrkSrcX_'//TRIM(FILE_NAME)
         call PutFile(TMP_NAME,BreakSourceX)
       ENDIF
       IF(OUT_BrkdisY)THEN
        TMP_NAME = TRIM(FDIR)//'BrkSrcY_'//TRIM(FILE_NAME)
         call PutFile(TMP_NAME,BreakSourceY)
       ENDIF


      IF(OUT_Time)THEN
        TMP_NAME = TRIM(FDIR)//'time_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,ARRTIME)
      ENDIF


     


! sediment

! end sediment

! foam

! end foam

     IF(OUT_TMP)THEN
        TMP_NAME = TRIM(FDIR)//'tmp_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,tmp4preview)
     ENDIF

101   continue

END SUBROUTINE PREVIEW

!-------------------------------------------------------------------------------------
!
!    PREVIEW_MEAN is subroutine for print-out of mean field data
!
!  HISTORY:
!    03/22/2016  Fengyan Shi
!-------------------------------------------------------------------------------------
SUBROUTINE PREVIEW_MEAN
     USE GLOBAL
     IMPLICIT NONE
     REAL(SP),DIMENSION(Mloc,Nloc) :: tmpout 

     CHARACTER(LEN=80)::FILE_NAME=' '
     CHARACTER(LEN=80)::FDIR=' '
     CHARACTER(LEN=80)::TMP_NAME=' '

     FDIR=TRIM(RESULT_FOLDER)

     ICOUNT_MEAN=ICOUNT_MEAN+1


        if (myid.eq.0)then
        WRITE(3,102)'PRINTING MEAN FILE', icount_mean
        WRITE(*,102)'PRINTING MEAN FILE', icount_mean
        endif


102     FORMAT(A20,I6)

	   itmp1=mod(icount_mean/10000,10)
	   itmp2=mod(icount_mean/1000,10)
	   itmp3=mod(icount_mean/100,10)
	   itmp4=mod(icount_mean/10,10)
	   itmp5=mod(icount_mean,10)

        write(file_name(1:1),'(I1)')itmp1
        write(file_name(2:2),'(I1)')itmp2
        write(file_name(3:3),'(I1)')itmp3
        write(file_name(4:4),'(I1)')itmp4
    	write(file_name(5:5),'(I1)')itmp5  

        IF(OUT_Umean)THEN
          TMP_NAME = TRIM(FDIR)//'umean_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,Umean)
          tmpout = P_mean / Max(Depth+ETAmean,MinDepthFrc)
          TMP_NAME = TRIM(FDIR)//'ulagm_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,tmpout)
        ENDIF
        IF(OUT_Vmean)THEN
          TMP_NAME = TRIM(FDIR)//'vmean_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,Vmean)
          tmpout = Q_mean / Max(Depth+ETAmean,MinDepthFrc)
          TMP_NAME = TRIM(FDIR)//'vlagm_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,tmpout)
        ENDIF
        IF(OUT_ETAmean)THEN
          TMP_NAME = TRIM(FDIR)//'etamean_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,ETAmean)
        ENDIF
        IF(OUT_WaveHeight)THEN
          TMP_NAME = TRIM(FDIR)//'Hrms_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,WaveHeightRMS)
          TMP_NAME = TRIM(FDIR)//'Havg_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,WaveHeightAve)
          TMP_NAME = TRIM(FDIR)//'Hsig_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,SigWaveHeight)
        ENDIF

 ! output FRCYsum,FRCYmean
 !        DxSxx,DySxy,DySyy,DxSxy,PgrdX,PgrdY,DxUUH,DyUVH,DyVVH,DxUVH

        IF(OUT_Radiation)THEN
          TMP_NAME = TRIM(FDIR)//'Sxx_'//TRIM(FILE_NAME)
          tmpout = UUmean-WWmean+0.5*9.8*ETA2mean
          call PutFile(TMP_NAME,tmpout)
          TMP_NAME = TRIM(FDIR)//'Sxy_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,UVmean)
          TMP_NAME = TRIM(FDIR)//'Syy_'//TRIM(FILE_NAME)
          tmpout = VVmean-WWmean+0.5*9.8*ETA2mean
          call PutFile(TMP_NAME,tmpout)

          TMP_NAME = TRIM(FDIR)//'DxSxx_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,DxSxx)
          TMP_NAME = TRIM(FDIR)//'DySxy_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,DySxy)
          TMP_NAME = TRIM(FDIR)//'DySyy_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,DySyy)
          TMP_NAME = TRIM(FDIR)//'DxSxy_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,DxSxy)
          TMP_NAME = TRIM(FDIR)//'PgrdX_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,PgrdX)
          TMP_NAME = TRIM(FDIR)//'PgrdY_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,PgrdY)
          TMP_NAME = TRIM(FDIR)//'DxUUH_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,DxUUH)
          TMP_NAME = TRIM(FDIR)//'DyUVH_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,DyUVH)
          TMP_NAME = TRIM(FDIR)//'DyVVH_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,DyVVH)
          TMP_NAME = TRIM(FDIR)//'DxUVH_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,DxUVH)
          TMP_NAME = TRIM(FDIR)//'FRCX_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,FRCXmean)
          TMP_NAME = TRIM(FDIR)//'FRCY_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,FRCYmean)
          TMP_NAME = TRIM(FDIR)//'BrkDissX_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,BreakDissX)
          TMP_NAME = TRIM(FDIR)//'BrkDissY_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,BreakDissY)

        ENDIF
                    

END SUBROUTINE PREVIEW_MEAN



!-------------------------------------------------------------------------------------
!
!    GetFile is subroutine for reading field data
!
!    HISTORY:
!    05/01/2010  Fengyan Shi
!    05/08/2017  Young-Kwang Choi
!-------------------------------------------------------------------------------------

SUBROUTINE GetFile(FILE,PHI)
     USE GLOBAL
     IMPLICIT NONE

     REAL(SP),DIMENSION(MGlob+2*Nghost,NGlob+2*Nghost) :: PHIGLOB
     CHARACTER(LEN=80) FILE
     REAL(SP),DIMENSION(Mloc,Nloc),INTENT(OUT) :: PHI

![-------ykchoi (08/May/2017)
     INTEGER :: irank, lenx, leny, lenxy, ireq
     INTEGER :: Nista, Niend, Njsta, Njend
     INTEGER :: istanum, iendnum, jstanum, jendnum
     INTEGER, ALLOCATABLE :: Nistas(:), Niends(:), Njstas(:), Njends(:)
     INTEGER :: istatus(mpi_status_size)
     REAL(SP), ALLOCATABLE :: xx(:,:)
! -------ykchoi (08/May/2017) ]

! TEMP

     if (myid.eq.0) then
        OPEN(1,FILE=TRIM(FILE))
        DO J=Nghost+1,NGlob+NGhost
           READ(1,*)(PHIGLOB(I,J),I=Nghost+1,MGlob+Nghost)
        ENDDO
        CLOSE(1)
! ghost cells
        DO I=Nghost+1,MGlob+Nghost
           DO J=1,Nghost
              PHIGLOB(I,J)=PHIGLOB(I,Nghost+1)
           ENDDO
           DO J=NGlob+Nghost+1,NGlob+2*Nghost
              PHIGLOB(I,J)=PHIGLOB(I,NGlob+Nghost)
           ENDDO
        ENDDO
        DO J=1,NGlob+2*Nghost
           DO I=1,Nghost
              PHIGLOB(I,J)=PHIGLOB(Nghost+1,J)
           ENDDO
           DO I=MGlob+Nghost+1,MGlob+2*Nghost
              PHIGLOB(I,J)=PHIGLOB(MGlob+Nghost,J)
           ENDDO
        ENDDO
     endif

![-------ykchoi (08/May/2017)
     Nista = iista + Nghost;
     Niend = iiend + Nghost;
     Njsta = jjsta + Nghost;
     Njend = jjend + Nghost;

     allocate( Nistas(nprocs), Niends(nprocs), Njstas(nprocs), Njends(nprocs) )

     call MPI_Gather( Nista, 1, MPI_INTEGER, Nistas, 1, MPI_INTEGER, &
                      0, MPI_COMM_WORLD, ier )
     call MPI_Gather( Niend, 1, MPI_INTEGER, Niends, 1, MPI_INTEGER, &
                      0, MPI_COMM_WORLD, ier )
     call MPI_Gather( Njsta, 1, MPI_INTEGER, Njstas, 1, MPI_INTEGER, &
                      0, MPI_COMM_WORLD, ier )
     call MPI_Gather( Njend, 1, MPI_INTEGER, Njends, 1, MPI_INTEGER, &
                      0, MPI_COMM_WORLD, ier )

     if( myid == 0 )then
	 PHI = PHIGLOB( 1:Mloc, 1:Nloc )
     endif

     do irank=1, px*py-1
	  if( myid == 0 ) then
	    istanum = Nistas(irank+1) - Nghost
	    iendnum = Niends(irank+1) + Nghost
	    jstanum = Njstas(irank+1) - Nghost
          jendnum = Njends(irank+1) + Nghost

	    lenx = iendnum - istanum + 1
	    leny = jendnum - jstanum + 1
	    lenxy = lenx*leny
	    allocate( xx(lenx, leny) )

	    xx = PHIGLOB( istanum:iendnum, jstanum:jendnum )
	    call mpi_isend( xx, lenxy, mpi_sp, irank, 1, mpi_comm_world, ireq, ier )
	    call mpi_wait( ireq, istatus, ier )
          deallocate( xx )

	  elseif( myid == irank ) then
	    
	    lenx = Niend-Nista+1+2*Nghost
	    leny = Njend-Njsta+1+2*Nghost
	    lenxy = lenx*leny

	    call mpi_irecv( PHI, lenxy, mpi_sp, 0, 1, mpi_comm_world, ireq, ier )
	    call mpi_wait( ireq, istatus, ier )

	  endif
     enddo

     deallocate( Nistas, Niends, Njstas, Njends )

! -------ykchoi (08/May/2017) ]

END SUBROUTINE Getfile




!-------------------------------------------------------------------------------------
!
!    PutFile is subroutine for print-out of field data
!
!    HISTORY:
!      05/01/2010  Fengyan Shi
!      05/06/2017  Young-Kwang Choi 
!-------------------------------------------------------------------------------------

SUBROUTINE PutFile(FILE_NAME,PHI)
     USE GLOBAL
     USE PARALLEL_FIELD_IO
     IMPLICIT NONE

     CHARACTER(LEN=80) FILE_NAME
     REAL(SP),DIMENSION(Mloc,Nloc),INTENT(IN) :: PHI

     CHARACTER(LEN=80)::TMP_NAME=' '

     SELECT CASE (TRIM(FIELD_IO_TYPE))
      CASE ('ASCII' , 'ascii')
         CALL PutFileASCII(FILE_NAME,PHI)
      CASE ('BINARY' , 'binary' )
         Call PutFileBinary(FILE_NAME,PHI)
      CASE DEFAULT
         !Defaults to ASCII case for non-valid input
         CALL PutFileASCII(FILE_NAME,PHI)
     END SELECT

END SUBROUTINE Putfile




! end vessel


! end sediment

! ----------

! end foam



