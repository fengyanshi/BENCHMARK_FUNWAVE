!------------------------------------------------------------------------------------
!
!      FILE mod_tide.F
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
!  TIDE_MODULE is a module to add tide/surge boundary conditions into wave simulation    
!
!  HISTORY :
!    11/26/2019 Fengyan Shi
!-------------------------------------------------------------------------------------

MODULE TIDE_MODULE
  USE PARAM
  USE GLOBAL,ONLY : Mloc,Nloc,Mloc1,Nloc1,Nghost,Ibeg,Iend,Jbeg,Jend,DX,DY, &
                    H,ETA,U,V,P,Q,MinDepth,MASK,DT,Gamma3,Depth,tmp4preview, &
                    ALPHA,BETA,MASK9,DepthX,DepthY,PERIODIC, &
                    UNDERTOW_U, UNDERTOW_V, ROLLER_SWITCH,ROLLER,Mglob,Nglob, &
                    TIME,Sponge_west_width,Sponge_east_width,Sponge_south_width, &
                    Sponge_north_width
  USE INPUT_READ

  USE GLOBAL,ONLY : myid,ier, npx,npy,PX,PY
  USE MPI

  IMPLICIT NONE
  SAVE

  CHARACTER(LEN=80) TideBcType
  CHARACTER(LEN=80) TideDataFileName
  CHARACTER(LEN=80) TideWestFileName,TideEastFileName,TideNorthFileName,TideSouthFileName
  CHARACTER(LEN=80) tmp_name
  LOGICAL :: TIDAL_BC_ABS=.FALSE.
  LOGICAL :: TIDAL_BC_GEN_ABS=.FALSE.
  LOGICAL :: TideEast=.TRUE.
  LOGICAL :: TideWest=.TRUE.
  LOGICAL :: TideNorth=.TRUE.
  LOGICAL :: TideSouth=.TRUE.
  REAL(SP) :: TideEast_U,TideEast_V,TideEast_ETA
  REAL(SP) :: TideWest_U,TideWest_V,TideWest_ETA
  REAL(SP) :: TideNorth_U,TideNorth_V,TideNorth_ETA
  REAL(SP) :: TideSouth_U,TideSouth_V,TideSouth_ETA
  REAL(SP),DIMENSION(:,:),ALLOCATABLE :: SPONGE_TIDE_WEST, &
              SPONGE_TIDE_EAST,SPONGE_TIDE_SOUTH,SPONGE_TIDE_NORTH
  INTEGER :: Iwidth,Ifile,Iwidth_tide
  REAL(SP) :: Time_tide_west_1,Eta_tide_west_1,U_tide_west_1,V_tide_west_1, &
              Time_tide_west_2,Eta_tide_west_2,U_tide_west_2,V_tide_west_2, &
              Time_tide_east_1,Eta_tide_east_1,U_tide_east_1,V_tide_east_1, &
              Time_tide_east_2,Eta_tide_east_2,U_tide_east_2,V_tide_east_2, &
              Time_tide_south_1,Eta_tide_south_1,U_tide_south_1,V_tide_south_1, &
              Time_tide_south_2,Eta_tide_south_2,U_tide_south_2,V_tide_south_2, &
              Time_tide_north_1,Eta_tide_north_1,U_tide_north_1,V_tide_north_1, &
              Time_tide_north_2,Eta_tide_north_2,U_tide_north_2,V_tide_north_2  


CONTAINS
  
SUBROUTINE TIDE_INITIAL
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,itmp5,SMALL,INPUT_FILE_NAME,WaveMaker

  USE GLOBAL,ONLY : iista,jjsta   

                    
  USE INPUT_READ
  IMPLICIT NONE

  CHARACTER(LEN=80)::FILE_NAME=' '
  CHARACTER(LEN=80)::TMP_NAME=' '
  INTEGER :: Ifile,ierr


! read  from input.txt
      FILE_NAME=INPUT_FILE_NAME

! width of sponge

      CALL READ_INTEGER(Iwidth_tide,FILE_NAME,'WaveMakerPointNum',ierr)

      IF(ierr==1)THEN
        Iwidth_tide = 30

      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'WaveMakerPointNum Default:  30 points'
         WRITE(3,'(A40)')'WaveMakerPointNum Default:  30 points'
      endif

       ENDIF     

 ! absorb only

      CALL READ_LOGICAL(TIDAL_BC_ABS,FILE_NAME,'TIDAL_BC_ABS',ierr)
      IF(ierr == 1)THEN
       TIDAL_BC_ABS = .FALSE. 

      if (myid.eq.0)then
       WRITE(3,'(A80)')'TIDAL_BC_ABS not defined, Default: False'
       WRITE(*,'(A80)')'TIDAL_BC_ABS not defined, Default: False'
      endif

      ENDIF

 ! generating and absorbing

      CALL READ_LOGICAL(TIDAL_BC_GEN_ABS,FILE_NAME,'TIDAL_BC_GEN_ABS',ierr)
      IF(ierr == 1)THEN
       TIDAL_BC_GEN_ABS = .FALSE. 
      ELSE
        IF(TIDAL_BC_GEN_ABS)THEN

        if (myid.eq.0)then
         WRITE(3,'(A80)')'TIDAL_BC_GEN_ABS is TRUE'
         WRITE(*,'(A80)')'TIDAL_BC_GEN_ABS is TRUE'
        endif

       ENDIF ! true
      ENDIF

!   tidal bc ---------
   IF(TIDAL_BC_GEN_ABS.OR.TIDAL_BC_ABS)THEN

! input tide details
      CALL READ_STRING(TideBcType,FILE_NAME,'TideBcType',ierr)
      IF(ierr==1)THEN
        TideBcType = 'CONSTANT'

      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'You dont specify TideBcType, use CONSTANT.'
         WRITE(3,'(A50)')'You dont specify TideBcType, use CONSTANT.'
      endif

      ENDIF

    IF(TideBcType(1:4)=='CONS')THEN


      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'You use CONSTANT TIDAL BC -----------------------'
         WRITE(3,'(A50)')'You use CONSTANT TIDAL BC -----------------------'
      endif


      CALL Tide_READ_CONSTANT
    ENDIF

    IF(TideBcType(1:4)=='DATA')THEN


      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'You use TIDAL BC DATA ---------------------------'
         WRITE(3,'(A50)')'You use TIDAL BC DATA ---------------------------'
      endif


      CALL Tide_READ_DATA_INIT
    ENDIF

      ALLOCATE(SPONGE_TIDE_WEST(Mloc,Nloc),SPONGE_TIDE_EAST(Mloc,Nloc), &
                SPONGE_TIDE_SOUTH(Mloc,Nloc),SPONGE_TIDE_NORTH(Mloc,Nloc))


      CALL TIDE_SPONGE

   ENDIF  ! end either tide abs or gen abs ---------------

    IF (TIDAL_BC_ABS) THEN


      if (myid.eq.0)then
       WRITE(3,'(A80)')'TIDAL_BC_ABS is defined'
       WRITE(*,'(A80)')'TIDAL_BC_ABS is defined'
      endif



      CALL REMOVE_SPONGE

    ENDIF ! end tide abs

END SUBROUTINE TIDE_INITIAL


SUBROUTINE TIDE_SPONGE
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,itmp5,SMALL

  USE GLOBAL,ONLY : iista,jjsta   

                    
  IMPLICIT NONE
  REAL(SP) :: ri,lim,Lstart,Lend,R_sp,A_sp
  
  Iwidth = Iwidth_tide ! specified in read or default = 30
  lim = 1.0_SP
  R_sp = 0.85_SP
  A_sp = 10.0_SP

! west

       do j = 1,Nloc
       do i = 1,Mloc

         ri = R_sp**(50*(i+npx*Mglob/px-1)/(Iwidth-1))

         SPONGE_TIDE_WEST(i,j) = max(A_sp**ri,lim)
       enddo
       enddo

! east
       do j = 1,Nloc
       do i = 1,Mloc

         ri = R_sp**(50*(Mloc-i+(px-npx-1)*Mglob/px)/(Iwidth-1))

         SPONGE_TIDE_EAST(i,j) = max(A_sp**ri,lim)
       enddo
       enddo

! south

       do j = 1,Nloc
       do i = 1,Mloc

         ri = R_sp**(50*(j+npy*Nglob/py-1)/(Iwidth-1))

         SPONGE_TIDE_SOUTH(i,j) = max(A_sp**ri,lim)
       enddo
       enddo

! north
       do j = 1,Nloc
       do i = 1,Mloc

         ri = R_sp**(50*(Nloc-j+(py-npy-1)*Nglob/py)/(Iwidth-1))

         SPONGE_TIDE_NORTH(i,j) = max(A_sp**ri,lim)
       enddo
       enddo

       SPONGE_TIDE_WEST = 1.0_SP/SPONGE_TIDE_WEST
       SPONGE_TIDE_EAST = 1.0_SP/SPONGE_TIDE_EAST
       SPONGE_TIDE_SOUTH = 1.0_SP/SPONGE_TIDE_SOUTH
       SPONGE_TIDE_NORTH = 1.0_SP/SPONGE_TIDE_NORTH



END SUBROUTINE TIDE_SPONGE

SUBROUTINE TIDE_BC
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,itmp5,SMALL

  USE GLOBAL,ONLY : iista,jjsta   

                    
  IMPLICIT NONE

   IF(TideWest)THEN
     DO J=1,Nloc
     DO I=1,Iwidth
       IF(MASK(I,J)==1)THEN
         ETA(I,J)=TideWest_ETA +(ETA(I,J)-TideWest_ETA)*SPONGE_TIDE_WEST(I,J)
         U(I,J)=TideWest_U +(U(I,J)-TideWest_U)*SPONGE_TIDE_WEST(I,J)
         V(I,J)=TideWest_V +(V(I,J)-TideWest_V)*SPONGE_TIDE_WEST(I,J)
       ENDIF
     ENDDO
     ENDDO
   ENDIF

   IF(TideEast)THEN
     DO J=1,Nloc
     DO I=Mloc-Iwidth+1,Mloc
       IF(MASK(I,J)==1)THEN
         ETA(I,J)=TideEast_ETA +(ETA(I,J)-TideEast_ETA)*SPONGE_TIDE_EAST(I,J)
         U(I,J)=TideEast_U +(U(I,J)-TideEast_U)*SPONGE_TIDE_EAST(I,J)
         V(I,J)=TideEast_V +(V(I,J)-TideEast_V)*SPONGE_TIDE_EAST(I,J)
       ENDIF
     ENDDO
     ENDDO
   ENDIF

   IF(TideSouth)THEN
     DO J=1,Iwidth
     DO I=1,Mloc
       IF(MASK(I,J)==1)THEN
         ETA(I,J)=TideSouth_ETA +(ETA(I,J)-TideSouth_ETA)*SPONGE_TIDE_SOUTH(I,J)
         U(I,J)=TideSouth_U +(U(I,J)-TideSouth_U)*SPONGE_TIDE_SOUTH(I,J)
         V(I,J)=TideSouth_V +(V(I,J)-TideSouth_V)*SPONGE_TIDE_SOUTH(I,J)
       ENDIF
     ENDDO
     ENDDO
   ENDIF

   IF(TideNorth)THEN
     DO J=Nloc-Iwidth+1,Nloc
     DO I=1,Mloc
       IF(MASK(I,J)==1)THEN
         ETA(I,J)=TideNorth_ETA +(ETA(I,J)-TideNorth_ETA)*SPONGE_TIDE_NORTH(I,J)
         U(I,J)=TideNorth_U +(U(I,J)-TideNorth_U)*SPONGE_TIDE_NORTH(I,J)
         V(I,J)=TideNorth_V +(V(I,J)-TideNorth_V)*SPONGE_TIDE_NORTH(I,J)
       ENDIF
     ENDDO
     ENDDO
   ENDIF

END SUBROUTINE TIDE_BC

SUBROUTINE Tide_READ_CONSTANT
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,itmp5,SMALL,INPUT_FILE_NAME,WaveMaker

  USE GLOBAL,ONLY : iista,jjsta   

                    
  USE INPUT_READ
  IMPLICIT NONE

  CHARACTER(LEN=80)::FILE_NAME=' '
  CHARACTER(LEN=80)::TMP_NAME=' '
  INTEGER :: Ifile,ierr
! read  from input.txt
      FILE_NAME=INPUT_FILE_NAME

!  west ---

      CALL READ_Float(TideWest_ETA,FILE_NAME,'TideWest_ETA',ierr)

      IF(ierr==1)THEN
        TideWest = .FALSE.

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideWest_ETA, use default: FALSE'
         WRITE(3,'(A80)')'You dont specify TideWest_ETA, use default: FALSE'
      endif

      ENDIF

      CALL READ_Float(TideWest_U,FILE_NAME,'TideWest_U',ierr)

      IF(ierr==1)THEN
        TideWest_U = ZERO

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideWest_U, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideWest_U, use default: 0'
      endif

      ENDIF

      CALL READ_Float(TideWest_V,FILE_NAME,'TideWest_V',ierr)

      IF(ierr==1)THEN
        TideWest_V = ZERO

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideWest_V, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideWest_V, use default: 0'
      endif

      ENDIF


!  east ---

      CALL READ_Float(TideEast_ETA,FILE_NAME,'TideEast_ETA',ierr)

      IF(ierr==1)THEN
        TideEast = .FALSE.

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideEast_ETA, use default: FALSE'
         WRITE(3,'(A80)')'You dont specify TideEast_ETA, use default: FALSE'
      endif

      ENDIF

      CALL READ_Float(TideEast_U,FILE_NAME,'TideEast_U',ierr)

      IF(ierr==1)THEN
        TideEast_U = ZERO

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideEast_U, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideEast_U, use default: 0'
      endif

      ENDIF

      CALL READ_Float(TideEast_V,FILE_NAME,'TideEast_V',ierr)

      IF(ierr==1)THEN
        TideEast_V = ZERO

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideEast_V, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideEast_V, use default: 0'
      endif

      ENDIF


!  south ---

      CALL READ_Float(TideSouth_ETA,FILE_NAME,'TideSouth_ETA',ierr)

      IF(ierr==1)THEN
        TideSouth = .FALSE.

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideSouth_ETA, use default: FALSE'
         WRITE(3,'(A80)')'You dont specify TideSouth_ETA, use default: FALSE'
      endif

      ENDIF

      CALL READ_Float(TideSouth_U,FILE_NAME,'TideSouth_U',ierr)

      IF(ierr==1)THEN
        TideSouth_U = ZERO

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideSouth_U, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideSouth_U, use default: 0'
      endif

      ENDIF

      CALL READ_Float(TideSouth_V,FILE_NAME,'TideSouth_V',ierr)

      IF(ierr==1)THEN
        TideSouth_V = ZERO

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideSouth_V, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideSouth_V, use default: 0'
      endif

      ENDIF

!  north ---

      CALL READ_Float(TideNorth_ETA,FILE_NAME,'TideNorth_ETA',ierr)

      IF(ierr==1)THEN
        TideNorth = .FALSE.

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideNorth_ETA, use default: FALSE'
         WRITE(3,'(A80)')'You dont specify TideNorth_ETA, use default: FALSE'
      endif

      ENDIF

      CALL READ_Float(TideNorth_U,FILE_NAME,'TideNorth_U',ierr)

      IF(ierr==1)THEN
        TideNorth_U = ZERO

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideNorth_U, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideNorth_U, use default: 0'
      endif

      ENDIF

      CALL READ_Float(TideNorth_V,FILE_NAME,'TideNorth_V',ierr)

      IF(ierr==1)THEN
        TideNorth_V = ZERO

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideNorth_V, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideNorth_V, use default: 0'
      endif

      ENDIF


END SUBROUTINE Tide_READ_CONSTANT


SUBROUTINE REMOVE_SPONGE
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,itmp5,SMALL,INPUT_FILE_NAME,WaveMaker

  USE GLOBAL,ONLY : iista,jjsta   

 
  IMPLICIT NONE

      IF(TideWest)THEN

        Sponge_west_width = ZERO

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You set TideWest BC,  Sponge_west is removed '
         WRITE(3,'(A80)')'You set TideWest BC,  Sponge_west is removed '
      endif

      ENDIF

      IF(TideEast)THEN

        Sponge_east_width = ZERO

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You set TideEast BC,  Sponge_east is removed '
         WRITE(3,'(A80)')'You set TideEast BC,  Sponge_east is removed '
      endif

      ENDIF

      IF(TideSouth)THEN

        Sponge_south_width = ZERO

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You set TideSouth BC,  Sponge_south is removed '
         WRITE(3,'(A80)')'You set TideSouth BC,  Sponge_south is removed '
      endif

      ENDIF

      IF(TideNorth)THEN

        Sponge_North_width = ZERO

      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You set TideNorth BC,  Sponge_North is removed '
         WRITE(3,'(A80)')'You set TideNorth BC,  Sponge_North is removed '
      endif

      ENDIF


      IF((.NOT.TideWest).AND.(.NOT.TideEast)  &
                        .AND.(.NOT.TideEast).AND.(.NOT.TideEast)) THEN
        TIDAL_BC_ABS = .FALSE.


      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You specify Tide but no BCs, set tide as false'
         WRITE(3,'(A80)')'You dont specify TideNorth_V, use default: 0'
      endif


      ENDIF ! end all false

END SUBROUTINE REMOVE_SPONGE

SUBROUTINE Tide_READ_DATA_INIT
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,itmp5,SMALL,INPUT_FILE_NAME,WaveMaker

  USE GLOBAL,ONLY : iista,jjsta   

                    
  USE INPUT_READ
  IMPLICIT NONE

  CHARACTER(LEN=80)::FILE_NAME=' '
  CHARACTER(LEN=80)::TMP_NAME=' '
  INTEGER :: Ifile,ierr
! read  from input.txt
      FILE_NAME=INPUT_FILE_NAME


! west -----------
       CALL READ_STRING(TideWestFileName,FILE_NAME,'TideWestFileName',ierr)

! intel compiler does not work well with inquire, so change to ierr==1
        IF(ierr==1)THEN

        IF(MYID==0)  &
         WRITE(*,*) TRIM(TMP_NAME), ' No west tidal bc'

        TideWest = .False.
        ELSE
          TMP_NAME = TRIM(TideWestFileName)
        ENDIF
      
      IF(TideWest) THEN
        Ifile = 201
        OPEN(Ifile,FILE=TRIM(TMP_NAME))   !!! here
        READ(Ifile,'(A80)')  tmp_name 
        READ(Ifile,*) Time_tide_west_2, Eta_tide_west_2, U_tide_west_2, V_tide_west_2

        Time_tide_west_1 = Time_tide_west_2 
        Eta_tide_west_1 = Eta_tide_west_2 
        U_tide_west_1 = U_tide_west_2 
        V_tide_west_1 = V_tide_west_2


   IF(MYID==0)THEN
     WRITE(*,*)'Read tide west ------'
     WRITE(*,*)'Time,Eta,U,V = ', Time_tide_west_2, Eta_tide_west_2, U_tide_west_2, V_tide_west_2
     WRITE(3,*)'Read tide ------'
     WRITE(3,*)'Time,Eta,U,V = ', Time_tide_west_2, Eta_tide_west_2, U_tide_west_2, V_tide_west_2
   ENDIF


     ENDIF ! tidewest


! east -----------
       CALL READ_STRING(TideEastFileName,FILE_NAME,'TideEastFileName',ierr)
        IF(ierr==1)THEN

        IF(MYID==0)  &
         WRITE(*,*) TRIM(TMP_NAME), ' No east tidal bc'

        TideEast = .False.
        ELSE
          TMP_NAME = TRIM(TideEastFileName)
        ENDIF
      
      IF(TideEast) THEN
        Ifile = 202
        OPEN(Ifile,FILE=TRIM(TMP_NAME))   !!! here
        READ(Ifile,'(A80)')  tmp_name 
        READ(Ifile,*) Time_tide_east_2, Eta_tide_east_2, U_tide_east_2, V_tide_east_2

        Time_tide_east_1 = Time_tide_east_2 
        Eta_tide_east_1 = Eta_tide_east_2 
        U_tide_east_1 = U_tide_east_2 
        V_tide_east_1 = V_tide_east_2


   IF(MYID==0)THEN
     WRITE(*,*)'Read tide east ------'
     WRITE(*,*)'Time,Eta,U,V = ', Time_tide_east_2, Eta_tide_east_2, U_tide_east_2, V_tide_east_2
     WRITE(3,*)'Read tide ------'
     WRITE(3,*)'Time,Eta,U,V = ', Time_tide_east_2, Eta_tide_east_2, U_tide_east_2, V_tide_east_2
   ENDIF


     ENDIF ! tideeast

! south -----------
       CALL READ_STRING(TideSouthFileName,FILE_NAME,'TideSouthFileName',ierr)
        IF(ierr==1)THEN

        IF(MYID==0)  &
         WRITE(*,*) TRIM(TMP_NAME), ' No south tidal bc'

        TideSouth = .False.
        ELSE
          TMP_NAME = TRIM(TideSouthFileName)
        ENDIF
      
      IF(TideSouth) THEN
        Ifile = 203
        OPEN(Ifile,FILE=TRIM(TMP_NAME))   !!! here
        READ(Ifile,'(A80)')  tmp_name 
        READ(Ifile,*) Time_tide_south_2, Eta_tide_south_2, U_tide_south_2, V_tide_south_2

        Time_tide_south_1 = Time_tide_south_2 
        Eta_tide_south_1 = Eta_tide_south_2 
        U_tide_south_1 = U_tide_south_2 
        V_tide_south_1 = V_tide_south_2


   IF(MYID==0)THEN
     WRITE(*,*)'Read tide south ------'
     WRITE(*,*)'Time,Eta,U,V = ', Time_tide_south_2, Eta_tide_south_2, U_tide_south_2, V_tide_south_2
     WRITE(3,*)'Read tide ------'
     WRITE(3,*)'Time,Eta,U,V = ', Time_tide_south_2, Eta_tide_south_2, U_tide_south_2, V_tide_south_2
   ENDIF


     ENDIF ! tidesouth

! north -----------
       CALL READ_STRING(TideNorthFileName,FILE_NAME,'TideNorthFileName',ierr)
        IF(ierr==1)THEN

        IF(MYID==0)  &
         WRITE(*,*) TRIM(TMP_NAME), ' No north tidal bc'

        TideNorth = .False.
        ELSE
          TMP_NAME = TRIM(TideNorthFileName)
        ENDIF
      
      IF(TideNorth) THEN
        Ifile = 204
        OPEN(Ifile,FILE=TRIM(TMP_NAME))   !!! here
        READ(Ifile,'(A80)')  tmp_name 
        READ(Ifile,*) Time_tide_north_2, Eta_tide_north_2, U_tide_north_2, V_tide_north_2

        Time_tide_north_1 = Time_tide_north_2 
        Eta_tide_north_1 = Eta_tide_north_2 
        U_tide_north_1 = U_tide_north_2 
        V_tide_north_1 = V_tide_north_2


   IF(MYID==0)THEN
     WRITE(*,*)'Read tide north ------'
     WRITE(*,*)'Time,Eta,U,V = ', Time_tide_north_2, Eta_tide_north_2, U_tide_north_2, V_tide_north_2
     WRITE(3,*)'Read tide ------'
     WRITE(3,*)'Time,Eta,U,V = ', Time_tide_north_2, Eta_tide_north_2, U_tide_north_2, V_tide_north_2
   ENDIF


     ENDIF ! tidenorth


END SUBROUTINE Tide_READ_DATA_INIT

SUBROUTINE TIDE_DATA
    USE GLOBAL,ONLY : tmp1,tmp2,SMALL,TIME,ZERO,DT
        INTEGER :: Ifile


! west ------------------
      IF(TideWest) THEN
        Ifile = 201  ! west

        IF(TIME>Time_tide_west_1.AND.TIME>Time_tide_west_2) THEN

        Time_tide_west_1 = Time_tide_west_2 
        Eta_tide_west_1 = Eta_tide_west_2 
        U_tide_west_1 = U_tide_west_2 
        V_tide_west_1 = V_tide_west_2

        DO WHILE (Time_tide_west_2.LT.TIME+DT)
          READ(Ifile,*,END=120)  Time_tide_west_2, Eta_tide_west_2, U_tide_west_2, V_tide_west_2
        END DO


   IF(MYID==0)THEN
     WRITE(*,*)'Read west ------'
     WRITE(*,*)'Time,Eta,U,V = ', Time_tide_west_2, Eta_tide_west_2, U_tide_west_2, V_tide_west_2
     WRITE(3,*)'Read tide ------'
     WRITE(3,*)'Time,Eta,U,V = ', Time_tide_west_2, Eta_tide_west_2, U_tide_west_2, V_tide_west_2
   ENDIF


     ENDIF ! end read time

! interpolate
    tmp2=ZERO
    tmp1=ZERO

    IF(TIME>Time_tide_west_1)THEN
      IF(Time_tide_west_1.EQ.Time_tide_west_2)THEN
        ! no more data
        tmp2=ZERO
        tmp1=ZERO
      ELSE
      tmp2=(Time_tide_west_2-TIME) &
            /MAX(SMALL, ABS(Time_tide_west_2-Time_tide_west_1))
      tmp1=1.0_SP - tmp2
      ENDIF  ! no more data?
    ENDIF ! time>time_1


    TideWEST_U = U_tide_west_2*tmp1 +U_tide_west_1*tmp2
    TideWEST_V = V_tide_west_2*tmp1 +V_tide_west_1*tmp2
    TideWEST_ETA = ETA_tide_west_2*tmp1 +ETA_tide_west_1*tmp2

      ENDIF  ! west

! east ------------------
      IF(TideEast) THEN
        Ifile = 202  ! east

        IF(TIME>Time_tide_east_1.AND.TIME>Time_tide_east_2) THEN

        Time_tide_east_1 = Time_tide_east_2 
        Eta_tide_east_1 = Eta_tide_east_2 
        U_tide_east_1 = U_tide_east_2 
        V_tide_east_1 = V_tide_east_2

        DO WHILE (Time_tide_east_2.LT.TIME+DT)
          READ(Ifile,*,END=120)  Time_tide_east_2, Eta_tide_east_2, U_tide_east_2, V_tide_east_2
        END DO


   IF(MYID==0)THEN
     WRITE(*,*)'Read east ------'
     WRITE(*,*)'Time,Eta,U,V = ', Time_tide_east_2, Eta_tide_east_2, U_tide_east_2, V_tide_east_2
     WRITE(3,*)'Read tide ------'
     WRITE(3,*)'Time,Eta,U,V = ', Time_tide_east_2, Eta_tide_east_2, U_tide_east_2, V_tide_east_2
   ENDIF


     ENDIF ! end read time

! interpolate
    tmp2=ZERO
    tmp1=ZERO

    IF(TIME>Time_tide_east_1)THEN
      IF(Time_tide_east_1.EQ.Time_tide_east_2)THEN
        ! no more data
        tmp2=ZERO
        tmp1=ZERO
      ELSE
      tmp2=(Time_tide_east_2-TIME) &
            /MAX(SMALL, ABS(Time_tide_east_2-Time_tide_east_1))
      tmp1=1.0_SP - tmp2
      ENDIF  ! no more data?
    ENDIF ! time>time_1

    TideEAST_U = U_tide_east_2*tmp1 +U_tide_east_1*tmp2
    TideEAST_V = V_tide_east_2*tmp1 +V_tide_east_1*tmp2
    TideEAST_ETA = ETA_tide_east_2*tmp1 +ETA_tide_east_1*tmp2

      ENDIF  ! east

! south ------------------
      IF(TideSouth) THEN
        Ifile = 203  ! east

        IF(TIME>Time_tide_south_1.AND.TIME>Time_tide_south_2) THEN

        Time_tide_south_1 = Time_tide_south_2 
        Eta_tide_south_1 = Eta_tide_south_2 
        U_tide_south_1 = U_tide_south_2 
        V_tide_south_1 = V_tide_south_2

        DO WHILE (Time_tide_south_2.LT.TIME+DT)
          READ(Ifile,*,END=120)  Time_tide_south_2, Eta_tide_south_2, U_tide_south_2, V_tide_south_2
        END DO


   IF(MYID==0)THEN
     WRITE(*,*)'Read south ------'
     WRITE(*,*)'Time,Eta,U,V = ', Time_tide_south_2, Eta_tide_south_2, U_tide_south_2, V_tide_south_2
     WRITE(3,*)'Read tide ------'
     WRITE(3,*)'Time,Eta,U,V = ', Time_tide_south_2, Eta_tide_south_2, U_tide_south_2, V_tide_south_2
   ENDIF


     ENDIF ! end read time

! interpolate
    tmp2=ZERO
    tmp1=ZERO

    IF(TIME>Time_tide_south_1)THEN
      IF(Time_tide_south_1.EQ.Time_tide_south_2)THEN
        ! no more data
        tmp2=ZERO
        tmp1=ZERO
      ELSE
      tmp2=(Time_tide_south_2-TIME) &
            /MAX(SMALL, ABS(Time_tide_south_2-Time_tide_south_1))
      tmp1=1.0_SP - tmp2
      ENDIF  ! no more data?
    ENDIF ! time>time_1

    TideSOUTH_U = U_tide_south_2*tmp1 +U_tide_south_1*tmp2
    TideSOUTH_V = V_tide_south_2*tmp1 +V_tide_south_1*tmp2
    TideSOUTH_ETA = ETA_tide_south_2*tmp1 +ETA_tide_south_1*tmp2

      ENDIF  ! south

! north ------------------
      IF(TideNorth) THEN
        Ifile = 204  ! north

        IF(TIME>Time_tide_north_1.AND.TIME>Time_tide_north_2) THEN

        Time_tide_north_1 = Time_tide_north_2 
        Eta_tide_north_1 = Eta_tide_north_2 
        U_tide_north_1 = U_tide_north_2 
        V_tide_north_1 = V_tide_north_2

        DO WHILE (Time_tide_north_2.LT.TIME+DT)
          READ(Ifile,*,END=120)  Time_tide_north_2, Eta_tide_north_2, U_tide_north_2, V_tide_north_2
        END DO


   IF(MYID==0)THEN
     WRITE(*,*)'Read north ------'
     WRITE(*,*)'Time,Eta,U,V = ', Time_tide_north_2, Eta_tide_north_2, U_tide_north_2, V_tide_north_2
     WRITE(3,*)'Read tide ------'
     WRITE(3,*)'Time,Eta,U,V = ', Time_tide_north_2, Eta_tide_north_2, U_tide_north_2, V_tide_north_2
   ENDIF


     ENDIF ! end read time

! interpolate
    tmp2=ZERO
    tmp1=ZERO

    IF(TIME>Time_tide_north_1)THEN
      IF(Time_tide_north_1.EQ.Time_tide_north_2)THEN
        ! no more data
        tmp2=ZERO
        tmp1=ZERO
      ELSE
      tmp2=(Time_tide_north_2-TIME) &
            /MAX(SMALL, ABS(Time_tide_north_2-Time_tide_north_1))
      tmp1=1.0_SP - tmp2
      ENDIF  ! no more data?
    ENDIF ! time>time_1

    TideNORTH_U = U_tide_north_2*tmp1 +U_tide_north_1*tmp2
    TideNORTH_V = V_tide_north_2*tmp1 +V_tide_north_1*tmp2
    TideNORTH_ETA = ETA_tide_north_2*tmp1 +ETA_tide_north_1*tmp2

      ENDIF  ! north


120 CONTINUE  ! no more data 

END SUBROUTINE TIDE_DATA

 
END MODULE TIDE_MODULE



    




