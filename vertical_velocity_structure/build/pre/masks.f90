!------------------------------------------------------------------------------------
!
!      FILE masks.F
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
!    UPDATE_MASK is subroutine to update mask
!       note that mask also be updated in fluxes subroutine
!    HISTORY: 
!       05/28/2010 Fengyan Shi
!-------------------------------------------------------------------------------------
SUBROUTINE UPDATE_MASK
     USE GLOBAL

     IMPLICIT NONE
     INTEGER,DIMENSION(:,:),ALLOCATABLE :: MASKtmp,MASK9tmp
     REAL(SP),DIMENSION(:,:),ALLOCATABLE :: ETAtmp

     REAL(SP)::myvar


     ALLOCATE(MASKtmp(Mloc,Nloc),MASK9tmp(Mloc,Nloc))
     ALLOCATE(ETAtmp(Mloc,Nloc))

! for the serial code, MASK at ghost cells keep no change


     call phi_int_exch(MASK)
     call phi_int_exch(MASK9)


     MASKtmp=MASK
! version 3.4 -------- 
!     ETAtmp=ETA
     Dmass = ZERO



     DO J=Jbeg-2,Jend+2
     DO I=Ibeg-2,Iend+2
! flood
     IF(MASK_STRUC(I,J)==1)THEN
       IF(MASK(I,J)<1)THEN
         ! left

           ! this is for rainfall case $$$
         IF(Eta(I,J)>-Depth(I,J))THEN
          MASKtmp(I,J)=1
         ENDIF

         IF(MASK(I-1,J)==1.AND.Eta(I-1,J)>Eta(I,J))THEN
           MASKtmp(I,J)=1
         ENDIF
         ! right
         IF(MASK(I+1,J)==1.AND.Eta(I+1,J)>Eta(I,J))THEN
           MASKtmp(I,J)=1
         ENDIF
         ! bottom
         IF(MASK(I,J-1)==1.AND.Eta(I,J-1)>Eta(I,j))THEN
           MASKtmp(I,J)=1
         ENDIF
         ! top
         IF(MASK(I,J+1)==1.AND.Eta(I,J+1)>Eta(I,j))THEN
           MASKtmp(I,J)=1
         ENDIF

! drying
       ELSE
         IF(Eta(I,J)<-Depth(I,J))THEN
          MASKtmp(I,J)=0
! version 3.4 -------- 
!          ETAtmp(I,J)=MinDepth-Depth(I,J)

          IF(I.GE.Ibeg.AND.I.LE.Iend.AND.J.GE.Jbeg.AND.J.LE.Jend)THEN

          Dmass = Dmass + (ETAtmp(I,J)-ETA(I,J))*DX*DY

          ENDIF ! ibeg-iend, jbeg,jend
         ENDIF    
       ENDIF
      ENDIF

! to avoid extreme depth gradient caused by depthx and depthy which were not
! treated in initialization, I reset depthx and depthy when drying 
! 01/21/2012
! HOWEVER, this truncation can affect the model accuracy as pointed by 
! Choi (private communication, in Thacker bowl test case, 07/03/2016). 
! In the following, I keep an option to use more accurate solution. 
! In Makefile, define -DIGNORE_BIG_SLOPE to ignore big slopes

        IF(MASK(I,J)<1)THEN
         DepthX(I,J)=Depth(I-1,J)
         DepthX(I+1,J)=Depth(I+1,J)
         DepthY(I,J)=Depth(I,J-1)
         DepthY(I,J+1)=Depth(I,J+1)
        ENDIF  


     ENDDO
     ENDDO


     MASK=MASKtmp
! version 3.4 -------- 
!     ETA=ETAtmp



     DO J=Jbeg-1,Jend+1
     DO I=Ibeg-1,Iend+1
      IF(VISCOSITY_BREAKING)THEN
        ! dont use mask9 for viscosity breaking
        MASK9tmp(I,J)= 1
      ELSE
        MASK9tmp(I,J)=MASK(I,J)*MASK(I-1,J)*MASK(I+1,J)  &
                *MASK(I+1,J+1)*MASK(I,J+1)*MASK(I-1,J+1) &
                *MASK(I+1,J-1)*MASK(I,J-1)*MASK(I-1,J-1) 
        IF(ABS(Eta(I,J))/MAX(DEPTH(I,J),MinDepthFrc)>SWE_ETA_DEP)THEN
         MASK9tmp(I,J)=ZERO
        ENDIF

       ENDIF ! end viscosity breaking

     ENDDO
     ENDDO

     MASK9 = MASK9tmp


  

     CALL PHI_INT_EXCH(MASK)
     CALL PHI_INT_EXCH(MASK9)


     DEALLOCATE(MASKtmp,MASK9tmp)
     DEALLOCATE(ETAtmp)


! end mass conservation

END SUBROUTINE UPDATE_MASK


