!------------------------------------------------------------------------------------
!
!      FILE mod_time_spectra.F
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
!  mod_time_spectra.F is a wavemaker module for time-dependent spectra   
!
!  HISTORY :
!    01/10/2024  Fengyan Shi
!
!  subroutines involved
!      io.F, main.F, mod_time_spectra.F (this one), mod_tide.F, sponge.F
!  model structure
!      * mainF, CALL READ_INPUT
!      * io.F, read WaveMaker = 'TIME_SPECTRA'
!      * main.F, CALL TIME_SPECTRA_INITIAL (mod_time_spectra.F)
!      * main.F, CALL TIDE_INITIAL (mod_tide.F), 
!                CALL TIME_SPECTRA_PROTECTION(mod_time_spectra.F)
!      * mod_time_spectra.F,  read spectra file and read first spectra data
!      * main.F, CALL TIDE_INITIAL (mod_tide.F)
!      * mod_tide.F, read TIDAL_BC_GEN_ABS, CALL Tide_READ_DATA_INIT, CALL TIDE_SPONGE
!      * main.F, CALL TIME_SPECTRA_PROTECTION
!      * main.F, in time loop, CALL TIME_SPECTRA_INTERPOLATION (mod_time_spectra.F)
!      * mod_time_spectra.F, read spectra data and interpolation
!        update Dep_Ser, CALL CALCULATE_DATA2D_Cm_Sm_TIME (mod_time_spectra.F)
!      * main.F, CALL ABSORBING_GENERATING_BC
!
!-------------------------------------------------------------------------------------


MODULE TIME_SPECTRA_MODULE
  USE PARAM
  USE GLOBAL,ONLY : Mloc,Nloc,Nghost,Ibeg,Iend,Jbeg,Jend,Mglob,Nglob, &
                    WidthWaveMaker,R_sponge_wavemaker,A_sponge_wavemaker, &
                    WaveCompFile,NumFreq,NumDir,&
                    Per_Ser,Theta_Ser,Phase_LEFT, Amp_Ser, &
                    Phase_Ser,&
                    Dep_Ser,Depth, INPUT_FILE_NAME,MASK,I,J,DX,DY,ZERO,SP, &
                    Grav,PI,TIME,Wave_Number_Ser,Stokes_Drift_Ser,Segma_Ser, &
                    Cm_eta,Sm_eta,Cm_u,Sm_u,Cm_v,Sm_v,PERIODIC,Beta_ref
  USE TIDE_MODULE
                     
  USE INPUT_READ

  USE GLOBAL,ONLY : myid,ier, npx,npy,PX,PY,iista,jjsta
  USE MPI

  IMPLICIT NONE
  SAVE
       REAL(SP), DIMENSION(:,:),ALLOCATABLE :: AmpData1,AmpData2

    REAL(SP) :: TimeSpectra1,TimeSpectra2
    CHARACTER (LEN=80) :: NameSpectra1,NameSpectra2
    CHARACTER(LEN=80)::FILE_NAME=' '
    CHARACTER(LEN=80)::SPECTRA_FILE =' ' 
    INTEGER :: PHASE_DATA     
       



    REAL(SP) :: myvar


CONTAINS
  
! READ Spectra

SUBROUTINE TIME_SPECTRA_INITIAL
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,SMALL,LARGE,INPUT_FILE_NAME
                    

  USE INPUT_READ
  IMPLICIT NONE

  INTEGER :: Ifile,ierr
  CHARACTER(LEN=80):: SpectraName
  CHARACTER(LEN=80) :: WHAT

! read precipitation from input.txt
      FILE_NAME=INPUT_FILE_NAME


      if (myid.eq.0) WRITE(3,*)'                                         '
      if (myid.eq.0) WRITE(3,*)'-------------- Time-dependent spectra INFO ----------'


! -------------------------

! spectra file
      CALL READ_STRING(SPECTRA_FILE,FILE_NAME,'SPECTRA_FILE',ierr)

      IF(ierr==1)THEN

        IF(MYID==0)  &
       WRITE(*,*) 'SPECTRA_FILE CANNOT BE FOUND. STOP'
       CALL MPI_FINALIZE (ier)
       STOP


      ELSE


      if (myid.eq.0) WRITE(3,'(A15,A50)')'SPECTRA_FILE:', SPECTRA_FILE


      ENDIF

! open file
  Ifile=400
  OPEN(Ifile,FILE=TRIM(SPECTRA_FILE))

! read file
         READ(Ifile,'(A80)')  WHAT ! title

      if (myid.eq.0) WRITE(*,*) WHAT

         READ(Ifile,'(A80)')  WHAT ! number of freq  and direction bins

      if (myid.eq.0) WRITE(*,*) WHAT

         READ(Ifile,*)  NumFreq,NumDir
           ALLOCATE (Amp_Ser(NumFreq,NumDir), AmpData1(NumFreq,NumDir),&
                     AmpData2(NumFreq,NumDir), &
                     Per_Ser(NumFreq),Theta_Ser(NumDir),Phase_LEFT(NumFreq,NumDir), &
                     Segma_Ser(NumFreq),Wave_Number_Ser(NumFreq))
           ALLOCATE(Cm_eta(Mloc,Nloc,Numfreq),Sm_eta(Mloc,Nloc,Numfreq), &
                    Cm_u(Mloc,Nloc,Numfreq),Sm_u(Mloc,Nloc,Numfreq),&
                    Cm_v(Mloc,Nloc,Numfreq),Sm_v(Mloc,Nloc,Numfreq) )

         READ(Ifile,'(A80)')  WHAT ! frequency bins

      if (myid.eq.0) WRITE(*,*) WHAT

       DO J=1,NumFreq
          READ(Ifile,*)Per_Ser(J)  ! read in as frequency
          Per_Ser(J)=1.0_SP/Per_Ser(J)
       ENDDO
         READ(Ifile,'(A80)')  WHAT ! direction bins

      if (myid.eq.0) WRITE(*,*) WHAT

       DO I=1,NumDir
          READ(Ifile,*)Theta_Ser(I)
       ENDDO
         READ(Ifile,'(A80)')  WHAT ! phases, the following can be empty 

      if (myid.eq.0) WRITE(*,*) WHAT

         READ(Ifile,*) Phase_Data
         IF (Phase_Data == 1) THEN
           DO I=1,NumDir
             READ(Ifile,*)(Phase_LEFT(J,I),J=1,NumFreq)
           ENDDO
           DO J=1,NumFreq
           DO I=1,NumDir
             Phase_LEFT(J,I)=Phase_LEFT(J,I)*3.1415926/180.0_SP
           ENDDO
           ENDDO

         ELSE
           DO J=1,NumFreq
           DO I=1,NumDir

              Phase_LEFT(J,I)=rand(0)*2.0_SP*3.1415926

           ENDDO
           ENDDO
         ENDIF ! phase true


         READ(Ifile,'(A80)')  WHAT ! t, file name

      if (myid.eq.0) WRITE(*,*) WHAT

         READ(Ifile,*)  TimeSpectra2
         READ(Ifile,'(A80)')  NameSpectra2

         TimeSpectra1 = TimeSpectra2
         NameSpectra1 = NameSpectra2


   IF(MYID==0)THEN
   WRITE(3,*) 'Initial Time, FileName: ', TimeSpectra2,TRIM(NameSpectra2)
   WRITE(*,*) 'Initial Time, FileName: ', TimeSpectra2,TRIM(NameSpectra2)
   ENDIF


! read data
  Ifile=401
  OPEN(Ifile,FILE=TRIM(NameSpectra2))

       DO I=1,NumDir
         READ(Ifile,*)(AmpData2(J,I),J=1,NumFreq)
!print*,I,AmpData2(J,I)
       ENDDO


  CLOSE(Ifile)


End SUBROUTINE TIME_SPECTRA_INITIAL

SUBROUTINE TIME_SPECTRA_PROTECTION
  USE GLOBAL,ONLY : tmp1,tmp2,SMALL,TIME,ZERO

  USE GLOBAL,ONLY : myid,ier, npx,npy,PX,PY
  USE MPI

  USE TIDE_MODULE
  
  IF (TIDAL_BC_GEN_ABS) THEN
    ! do nothing
  ELSE

   IF(MYID==0)THEN
   WRITE(3,*) 'You should specify TIDAL_BC_GEN_ABS, STOP!'
   WRITE(*,*) 'You should specify TIDAL_BC_GEN_ABS, STOP!'
   ENDIF
          call MPI_FINALIZE ( ier )

  STOP
  ENDIF

  IF (TideBcType(1:4)=='DATA') THEN
    ! do nothing
  ELSE

   IF(MYID==0)THEN
   WRITE(3,*) 'You should specify TideBcType = DATA, STOP!'
   WRITE(*,*) 'You should specify TideBcType, STOP!'
   ENDIF
          call MPI_FINALIZE ( ier )

  STOP
  ENDIF

END SUBROUTINE TIME_SPECTRA_PROTECTION

SUBROUTINE TIME_SPECTRA_INTERPOLATION
  USE GLOBAL,ONLY : tmp1,tmp2,SMALL,TIME,ZERO
  USE INPUT_READ
  IMPLICIT NONE
  INTEGER :: Ifile,ierr,I,J
  REAL(SP) :: rII,rJJ
  INTEGER :: IOstatus


    IF(TIME>TimeSpectra1.AND.TIME>TimeSpectra2) THEN

         TimeSpectra1=TimeSpectra2
         NameSpectra1=NameSpectra2
         AmpData1=AmpData2

    Ifile = 400

    READ(Ifile,*,IOSTAT=IOstatus)  TimeSpectra2

    IF(IOstatus< 0)GOTO 120

    READ(Ifile,*)  NameSpectra2


   IF(MYID==0)THEN
   WRITE(3,*) 'READ Spectra, Time, FileName: ', TimeSpectra2,TRIM(NameSpectra2)
   WRITE(*,*) 'READ Spectra, Time, FileName: ', TimeSpectra2,TRIM(NameSpectra2)
   ENDIF


! read data
  Ifile=401

  OPEN(Ifile,FILE=TRIM(NameSpectra2))
       DO I=1,NumDir
         READ(Ifile,*)(AmpData2(J,I),J=1,NumFreq)
!print*,I,Amp_Ser(J,I)
       ENDDO
  CLOSE(Ifile)

    ENDIF ! end time > timeSpectra2

! intermpolation
    tmp2=ZERO
    tmp1=ZERO

    IF(TIME>TimeSpectra1)THEN
      IF(TimeSpectra1.EQ.TimeSpectra2)THEN
        ! no more data
        tmp2=ZERO
        tmp1=ZERO
      ELSE
      tmp2=(TimeSpectra2-TIME) &
            /MAX(SMALL, ABS(TimeSpectra2-TimeSpectra1))
      tmp1=1.0_SP - tmp2;
      ENDIF  ! no more data?
    ENDIF ! time>time_1

    Amp_Ser = AmpData2*tmp1 +AmpData1*tmp2

120 CONTINUE  ! no more data for vessel Kves

   Dep_Ser = TideWest_ETA + Depth(1,1)

     CALL CALCULATE_DATA2D_Cm_Sm_TIME


! debug -------

!    if(myid==0)then
!    open(500,file='tmp1.txt')
!    elseif(myid==1)then
!    open(500,file='tmp2.txt')
!    elseif(myid==2)then
!    open(500,file='tmp3.txt')
!    elseif(myid==3)then
!    open(500,file='tmp4.txt')
!    endif
!     do j=1,nloc
!       write(500,*)(PrecRateModel(i,j),i=1,mloc)
!     enddo
!    close(500)
!    stop
!    CALL MPI_FINALIZE (ier)
! debug over --------

END SUBROUTINE TIME_SPECTRA_INTERPOLATION

SUBROUTINE CALCULATE_DATA2D_Cm_Sm_TIME

     IMPLICIT NONE
     INTEGER :: Iter,KK,K,KKK
     REAL(SP) :: Celerity,Wave_Length,Fk,Fkdif,Zlev,Fact,X_maker,Y_maker
!     real(SP),DIMENSION(Mloc,Nloc) :: Ein2D,Din2D
!     real(SP),DIMENSION(Mloc,Nloc) :: Uin2D,Vin2D
     REAL(SP) :: tmp_eta,tmp_u,tmp_v
     REAL(SP) :: Theta_Per

! wave_number

   DO I=1,NumFreq
     Segma_Ser(I) = 2.0*pi/Per_Ser(I)
     Celerity = sqrt(Grav*Dep_Ser)
     Wave_Length = Celerity*Per_Ser(I)
     Wave_Number_Ser(I) = 2.0*pi/Wave_Length
     
     Iter = 0
55   Fk = Grav*Wave_Number_Ser(I)*tanh(Wave_Number_Ser(I)*Dep_Ser)-Segma_Ser(I)**2
     if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 65
     Fkdif = Grav*Wave_Number_Ser(I)*Dep_Ser*(1.0-tanh(Wave_Number_Ser(I)*Dep_Ser)**2)+  &
        Grav*tanh(Wave_Number_ser(I)*Dep_Ser) 
     Wave_Number_Ser(I) = Wave_Number_Ser(I)-Fk/Fkdif
     Iter = Iter+1
     goto 55
65   continue
     Wave_Length = 2.0*pi/Wave_Number_Ser(I)

    ENDDO ! end NumCompSer  

     Stokes_Drift_Ser = ZERO
     Fact = 1.0 

! Cm and Sm
     Cm_eta = ZERO
     Sm_eta = ZERO
     Cm_u = ZERO
     Sm_u = ZERO
     Cm_v = ZERO
     Sm_v = ZERO

     Zlev = ABS(1.0_SP+Beta_ref)*Dep_Ser
!     Zlev = Dep_Ser

     DO KK=1,NumFreq

         DO KKK=1,NumDir
          IF(PERIODIC)THEN
            CALL CalcPeriodicTheta(Wave_Number_Ser(KK),Theta_Ser(KKK),Theta_Per)
          ELSE
            Theta_Per = Theta_Ser(KKK)
          ENDIF 

       DO J=1,Nloc
       DO I=1,Mloc

         X_maker=( I-Ibeg+(iista-1) )*DX
         Y_maker=( J-Jbeg+(jjsta-1) )*DY




          Cm_eta(I,J,KK)=Cm_eta(I,J,KK)+Amp_Ser(KK,KKK)*COS( &
                   Wave_Number_Ser(KK)*SIN(Theta_Per)*Y_maker &
                  +Wave_Number_Ser(KK)*COS(Theta_Per)*X_maker )
          Sm_eta(I,J,KK)=Sm_eta(I,J,KK)+Amp_Ser(KK,KKK)*SIN( &
                   Wave_Number_Ser(KK)*SIN(Theta_Per)*Y_maker &
                  +Wave_Number_Ser(KK)*COS(Theta_Per)*X_maker )
          Cm_u(I,J,KK)=Cm_u(I,J,KK)+Amp_Ser(KK,KKK) &
                       *Segma_Ser(KK)*cosh(Wave_Number_Ser(KK)*Zlev)  &
                       /sinh(Wave_Number_Ser(KK)*Dep_Ser)  &
                   *COS(Theta_Per)*COS( &
                   Wave_Number_Ser(KK)*SIN(Theta_Per)*Y_maker &
                  +Wave_Number_Ser(KK)*COS(Theta_Per)*X_maker )
          Cm_v(I,J,KK)=Cm_v(I,J,KK)+Amp_Ser(KK,KKK) &
                       *Segma_Ser(KK)*cosh(Wave_Number_Ser(KK)*Zlev)  &
                       /sinh(Wave_Number_Ser(KK)*Dep_Ser)  &
                   *SIN(Theta_Per)*COS( &
                   Wave_Number_Ser(KK)*SIN(Theta_Per)*Y_maker &
                  +Wave_Number_Ser(KK)*COS(Theta_Per)*X_maker )                   
          Sm_u(I,J,KK)=Sm_u(I,J,KK)+Amp_Ser(KK,KKK) &
                       *Segma_Ser(KK)*cosh(Wave_Number_Ser(KK)*Zlev)  &
                       /sinh(Wave_Number_Ser(KK)*Dep_Ser)  &
                   *COS(Theta_Per)*SIN( &
                   Wave_Number_Ser(KK)*SIN(Theta_Per)*Y_maker &
                  +Wave_Number_Ser(KK)*COS(Theta_Per)*X_maker )
          Sm_v(I,J,KK)=Sm_v(I,J,KK)+Amp_Ser(KK,KKK) &
                       *Segma_Ser(KK)*cosh(Wave_Number_Ser(KK)*Zlev)  &
                       /sinh(Wave_Number_Ser(KK)*Dep_Ser)  &
                   *SIN(Theta_Per)*SIN( &
                   Wave_Number_Ser(KK)*SIN(Theta_Per)*Y_maker &
                  +Wave_Number_Ser(KK)*COS(Theta_Per)*X_maker )

         ENDDO

       ENDDO
       ENDDO
      ENDDO


END SUBROUTINE CALCULATE_DATA2D_Cm_Sm_TIME


END MODULE TIME_SPECTRA_MODULE


