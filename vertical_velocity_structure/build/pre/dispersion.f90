!------------------------------------------------------------------------------------
!
!      FILE dispersion.F
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
!    CAL_DISPERSION is subroutine to calculation dispersion terms
!    so far V^4 and V^1
!    called by
!       MAIN
!    call DERIVATIVE_XX
!         DERIVATIVE_XY
!    
!    HISTORY: 
!      05/01/2010 Fengyan Shi
!      10/14/2012 Fengyan Shi, added coupling bc,
!                              change derivative_xx_high to second order
!                              according to Harris suggestion
!      08/06/2015 - 08/18/2015 Young-Kwang Choi, modified t-derivatives
!                   corrected U1p,V1p, V2, V3 and omega_1 terms
!
!-------------------------------------------------------------------------------------
SUBROUTINE CAL_DISPERSION
     USE GLOBAL
     IMPLICIT NONE

     REAL(SP),Dimension(Mloc,Nloc) :: DU,DV,DUt,DVt
     REAL(SP) :: UxxVxy,UxyVyy,HUxxHVxy,HUxyHVyy, &
                 UxxVxy_x,UxxVxy_y,UxyVyy_x,UxyVyy_y, &
                 HUxxHVxy_x,HUxxHVxy_y,HUxyHVyy_x,HUxyHVyy_y, &
                 rh,rhx,rhy,reta,ken1,ken2,ken3,ken4,ken5

     REAL(SP) :: omega_0,omega_1
     REAL(SP),Dimension(Mloc,Nloc) :: omega            

    
! uxx
    CALL DERIVATIVE_XX(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,U,Uxx)
! uxy
    CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DY,U,Uxy)
! vxy
    CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DY,V,Vxy)
! vyy
    CALL DERIVATIVE_YY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,V,Vyy)


! gamma2.ne.0
    IF(Gamma2>ZERO)THEN
     CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,U,Ux)
     CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,V,Vx)
     CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,U,Uy)
     CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,V,Vy)
     CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,Eta,ETAx)
     CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,Eta,ETAy)
    ELSEIF(SHOW_BREAKING)THEN

     CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,Eta,ETAx)
     CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,Eta,ETAy)
    ENDIF

! DU DV
     DO J=1,Nloc-1
     DO I=1,Mloc-1
       DU(I,J)=Max(Depth(I,J),MinDepthFrc)*U(I,J)
       DV(I,J)=Max(Depth(I,J),MinDepthFrc)*V(I,J)

! ykchoi (15. 08. 18)
! Computation of Etat for ( U1p )_t in conservative form of Shi et al. (2012)
! and for viscosity of wave maker
! 05/25/2018, ykchoi pointed a bug: J=1,Nloc, I=1,Mloc


       ETAT(I,J)=-(P(I+1,J)-P(I,J))/DX-(Q(I,J+1)-Q(I,J))/DY

     ENDDO
     ENDDO 

! ETAT

    IF(Gamma2>ZERO)THEN


! as pointed by Steve Brandt, Ut and Vt should be evaluated outside of RK loop.
! Because it is really small term, I temporarily keep this form but need more tests
! to see if or not affect results. It is also important for hot start

! ykchoi (15. 08. 06)
! Modification : Ut, Vt

       DO J=1,Nloc
       DO I=1,Mloc

	  Ut(I,J) = (U(I,J)-U0(I,J)) / DT   !ykchoi
	  Vt(I,J) = (V(I,J)-V0(I,J)) / DT   !ykchoi

        DUt(I,J)=Max(Depth(I,J),MinDepthFrc)*Ut(I,J)
        DVt(I,J)=Max(Depth(I,J),MinDepthFrc)*Vt(I,J)
       ENDDO
       ENDDO

    ELSEIF(SHOW_BREAKING .OR. WAVEMAKER_VIS)THEN


    ENDIF

! DUxx
    CALL DERIVATIVE_XX(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DU,DUxx)
! DUxy
    CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DY,DU,DUxy)
! DVxy
    CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DY,DV,DVxy)
! DVyy
    CALL DERIVATIVE_YY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,DV,DVyy)

      

    IF(Gamma2>ZERO)THEN
     CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DU,DUx)
     CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DV,DVx)
     CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,DU,DUy)
     CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,DV,DVy)
     CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,Ut,Utx)
     CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,Vt,Vty)

     CALL DERIVATIVE_XX(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,Ut,Utxx)
     CALL DERIVATIVE_YY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,Vt,Vtyy)
     CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DY,Ut,Utxy)
     CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DY,Vt,Vtxy)

     CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DUt,DUtx)
     CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,DVt,DVty)

     CALL DERIVATIVE_XX(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DUt,DUtxx)
     CALL DERIVATIVE_YY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,DVt,DVtyy)
     CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DY,DUt,DUtxy)
     CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DY,DVt,DVtxy)
    ENDIF


! this may affect parallel version
! I added coupling boundary 10/14/2012
!  add left_bc wavemaker 09/12/2017

!  boundary conditions

    if(n_west.eq.MPI_PROC_NULL) then

    IF (WaveMaker(1:11)=='LEFT_BC_IRR')THEN
     ! do nothing
    ELSE
     DO J=1,Nloc
       Uxy(Ibeg,J)=ZERO
       DUxy(Ibeg,J)=ZERO
       Vxy(Ibeg,J)=ZERO
       DVxy(Ibeg,J)=ZERO
       Utxy(Ibeg,J)=ZERO
       DUtxy(Ibeg,J)=ZERO
       Vtxy(Ibeg,J)=ZERO
       DVtxy(Ibeg,J)=ZERO
     ENDDO
    ENDIF ! left_bc wavemaker




    endif  



    if(n_east.eq.MPI_PROC_NULL) then

     DO J=1,Nloc
       Uxy(Iend,J)=ZERO
       DUxy(Iend,J)=ZERO
       Vxy(Iend,J)=ZERO
       DVxy(Iend,J)=ZERO
       Utxy(Iend,J)=ZERO
       DUtxy(Iend,J)=ZERO
       Vtxy(Iend,J)=ZERO
       DVtxy(Iend,J)=ZERO
     ENDDO 

    endif  

  

    if(n_suth.eq.MPI_PROC_NULL) then

     DO I=1,Mloc
       Uxy(I,Jbeg)=ZERO
       DUxy(I,Jbeg)=ZERO
       Vxy(I,Jbeg)=ZERO
       DVxy(I,Jbeg)=ZERO
       Utxy(I,Jbeg)=ZERO
       DUtxy(I,Jbeg)=ZERO
       Vtxy(I,Jbeg)=ZERO
       DVtxy(I,Jbeg)=ZERO
     ENDDO   

    endif  



    if(n_nrth.eq.MPI_PROC_NULL) then

     DO I=1,Mloc
       Uxy(I,Jend)=ZERO
       DUxy(I,Jend)=ZERO
       Vxy(I,Jend)=ZERO
       DVxy(I,Jend)=ZERO
       Utxy(I,Jend)=ZERO
       DUtxy(I,Jend)=ZERO
       Vtxy(I,Jend)=ZERO
       DVtxy(I,Jend)=ZERO
     ENDDO 

    endif  



    CALL EXCHANGE_DISPERSION

     
! calculate V1p  without nonlinear dispersion
     DO J=1,Nloc
     DO I=1,Mloc

       U4(I,J)=(1.0_SP/3.0_SP-Beta_1+0.5_SP*Beta_1*Beta_1)*DEPTH(I,J)*DEPTH(I,J)*(Uxx(I,J)+Vxy(I,J)) &
                +(Beta_1-1.0_SP/2.0_SP)*DEPTH(I,J)*(DUxx(I,J)+DVxy(I,J))
       V4(I,J)=(1.0_SP/3.0_SP-Beta_1+0.5_SP*Beta_1*Beta_1)*DEPTH(I,J)*DEPTH(I,J)*(Uxy(I,J)+Vyy(I,J)) &
                +(Beta_1-1.0_SP/2.0_SP)*DEPTH(I,J)*(DUxy(I,J)+DVyy(I,J))      

       !------[ykchoi (04/14/2017)
       IF(gamma2>ZERO)THEN
	   UxxVxy = Uxx(I,J) + Vxy(I,J)
	   UxyVyy = Uxy(I,J) + Vyy(I,J)

	   HUxxHVxy = DUxx(I,J) + DVxy(I,J)
	   HUxyHVyy = DUxy(I,J) + DVyy(I,J)

	   rh = Depth(I,J)
	   reta = Eta(I,J)

  	   ken1 = ( 1.0_SP/6.0_SP - Beta_1 + Beta_1*Beta_1 )*rh*reta*Beta_2   &
	         + ( 1.0_SP/2.0_SP*Beta_1*Beta_1 - 1.0_SP/6.0_SP )*reta*reta*Beta_2*Beta_2
	   ken2 = ( Beta_1 - 1.0_SP/2.0_SP )*reta*Beta_2

	   U4(I,J) = U4(I,J) + gamma2*MASK9(I,J)*( ken1*UxxVxy + ken2*HUxxHVxy )
	   V4(I,J) = V4(I,J) + gamma2*MASK9(I,J)*( ken1*UxyVyy + ken2*HUxyHVyy )
	 ENDIF
	 !------ykchoi (04/14/2017)]





! ykchoi( 15. 08. 06.)
! U1p, V1p terms are modified to 0.5_SP*(1.0_SP-Beta_1) --> 0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)
               

       !U1p(I,J)=0.5_SP*(1.0_SP-Beta_1)  & !ykchoi
	 U1p(I,J)=0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)  &
                *DEPTH(I,J)*DEPTH(I,J)  &
                *(Uxx(I,J)+Vxy(I,J)) &
               +(Beta_1-1.0_SP)*DEPTH(I,J)*(DUxx(I,J)+DVxy(I,J))
       
	 !V1p(I,J)=0.5_SP*(1.0_SP-Beta_1)  & !ykchoi
	 V1p(I,J)=0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)  &
                *DEPTH(I,J)*DEPTH(I,J)  &
                *(Uxy(I,J)+Vyy(I,J)) &
               +(Beta_1-1.0_SP)*DEPTH(I,J)*(DUxy(I,J)+DVyy(I,J))


     ENDDO
     ENDDO


       DO J=Jbeg,Jend
       DO I=Ibeg,Iend
        grdAx(I,J)=DUxx(I,J)+DVxy(I,J)
        grdAy(I,J)=DUxy(I,J)+DVyy(I,J)
        grdBx(I,J)=Uxx(I,J)+Vxy(I,J)
        grdBy(I,J)=Uxy(I,J)+Vyy(I,J)
       ENDDO
       ENDDO




     IF(gamma2>ZERO)THEN
       DO J=Jbeg,Jend
       DO I=Ibeg,Iend
! 
        UxxVxy=Uxx(I,J)+Vxy(I,J)
        UxyVyy=Uxy(I,J)+Vyy(I,J)
        UxxVxy_x=(Uxx(I+1,J)+Vxy(I+1,J)-Uxx(I-1,J)-Vxy(I-1,J))/2.0_SP/DX
        UxxVxy_y=(Uxx(I,J+1)+Vxy(I,J+1)-Uxx(I,J-1)-Vxy(I,J-1))/2.0_SP/DY
        UxyVyy_x=(Uxy(I+1,J)+Vyy(I+1,J)-Uxy(I-1,J)-Vyy(I-1,J))/2.0_SP/DX
        UxyVyy_y=(Uxy(I,J+1)+Vyy(I,J+1)-Uxy(I,J-1)-Vyy(I,J-1))/2.0_SP/DY

        HUxxHVxy=DUxx(I,J)+DVxy(I,J)
        HUxyHVyy=DUxy(I,J)+DVyy(I,J)
        HUxxHVxy_x=(DUxx(I+1,J)+DVxy(I+1,J)-DUxx(I-1,J)-DVxy(I-1,J))/2.0_SP/DX
        HUxxHVxy_y=(DUxx(I,J+1)+DVxy(I,J+1)-DUxx(I,J-1)-DVxy(I,J-1))/2.0_SP/DY
        HUxyHVyy_x=(DUxy(I+1,J)+DVyy(I+1,J)-DUxy(I-1,J)-DVyy(I-1,J))/2.0_SP/DX
        HUxyHVyy_y=(DUxy(I,J+1)+DVyy(I,J+1)-DUxy(I,J-1)-DVyy(I,J-1))/2.0_SP/DY

        rh=Depth(I,J)
        rhx=(Depth(I+1,J)-Depth(I-1,J))/2.0_SP/DX
        rhy=(Depth(I,J+1)-Depth(I,J-1))/2.0_SP/DY
        reta=Eta(I,J)

        U1pp(I,J)=-reta*Beta_2*ETAx(I,J)*Beta_2*(Utx(I,J)+Vty(I,J)) - 0.5_SP*reta*reta*Beta_2*Beta_2*(Utxx(I,J)+Vtxy(I,J))&
                  -ETAx(I,J)*Beta_2*(DUtx(I,J)+DVty(I,J)) -reta*Beta_2*(DUtxx(I,J)+DVtxy(I,J))

        V1pp(I,J)=-reta*Beta_2*ETAy(I,J)*Beta_2*(Utx(I,J)+Vty(I,J)) - 0.5_SP*reta*reta*Beta_2*Beta_2*(Utxy(I,J)+Vtyy(I,J))&
                  -ETAy(I,J)*Beta_2*(DUtx(I,J)+DVty(I,J)) -reta*Beta_2*(DUtxy(I,J)+DVtyy(I,J))
        
	  !ykchoi (2015. 08.06.)
	  !remaining terms in U1p, V1p in Shi et al. (2012) publication   
	  ken1 = Beta_1*( 1.0_SP - Beta_1 )*rh*ETAT(I,J)*Beta_2 - Beta_1*Beta_1*reta*Beta_2*ETAT(I,J)*Beta_2
	  ken2 = Beta_1*( 1.0_SP - Beta_1 )*rh*reta*Beta_2 - 0.5*Beta_1*Beta_1*reta*reta*Beta_2*Beta_2
	  ken3 = Beta_1*ETAT(I,J)*Beta_2
	  ken4 = Beta_1*reta*Beta_2

	  U1pp(I,J) = U1pp(I,J) - ken1*UxxVxy - ken2*( Utxx(I,J)+Vtxy(I,J) ) + ken3*HUxxHVxy + ken4*(DUtxx(I,J)+DVtxy(I,J))
	  V1pp(I,J) = V1pp(I,J) - ken1*UxyVyy - ken2*( Utxy(I,J)+Vtyy(I,J) ) + ken3*HUxyHVyy + ken4*(DUtxy(I,J)+DVtyy(I,J))

         ken1=(Beta_1-1.0_SP)*(rhx+ETAx(I,J))*Beta_2
         ken2=(Beta_1-1.0_SP)*(rh+reta)*Beta_2
         ken3=( (1.0_SP-Beta_1)*(1.0_SP-Beta_1)*rh*rhx*Beta_2*Beta_2-Beta_1*(1.0_SP-Beta_1)*(rhx*reta*Beta_2+rh*ETAx(I,J)*Beta_2) &
                    +(Beta_1*Beta_1-1.0_SP)*reta*ETAx(I,J)*Beta_2*Beta_2 )
         ken4=( 0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)*rh*rh*Beta_2*Beta_2-Beta_1*(1.0_SP-Beta_1)*rh*reta*Beta_2 &
                      +0.5_SP*(Beta_1*Beta_1-1.0_SP)*reta*reta*Beta_2*Beta_2 )
         ken5=( (1.0_SP-Beta_1)*(1.0_SP-Beta_1)*rh*rhy*Beta_2*Beta_2-Beta_1*(1.0_SP-Beta_1)*(rhy*reta*Beta_2+rh*ETAy(I,J)*Beta_2) &
                    +(Beta_1*Beta_1-1.0_SP)*reta*ETAy(I,J)*Beta_2*Beta_2 )

        U2(I,J)=ken1*(U(I,J)*HUxxHVxy+V(I,J)*HUxyHVyy) &
                +ken2*(Ux(I,J)*HUxxHVxy+U(I,J)*HUxxHVxy_x &
                    +Vx(I,J)*HUxyHVyy+V(I,J)*HUxyHVyy_x) &
                +ken3 & 
                   *(U(I,J)*UxxVxy+V(I,J)*UxyVyy) &
                +ken4  &
                   *(Ux(I,J)*UxxVxy+U(I,J)*UxxVxy_x+Vx(I,J)*UxyVyy+V(I,J)*UxyVyy_x) &
                +Beta_2*Beta_2*(DUx(I,J)+DVy(I,J)+reta*Beta_2*(Ux(I,J)+Vy(I,J)))  &
                   *(HUxxHVxy+ETAx(I,J)*Beta_2*(Ux(I,J)+Vy(I,J))+reta*Beta_2*UxxVxy)

        ! ykchoi (15. 08. 06.)
	  ! ken1 should be modified as following term.
        ken1=(Beta_1-1.0_SP)*(rhy+ETAy(I,J))*Beta_2

        V2(I,J)=ken1*(U(I,J)*HUxxHVxy+V(I,J)*HUxyHVyy) &
                +ken2*(Uy(I,J)*HUxxHVxy+U(I,J)*HUxxHVxy_y &
                    +Vy(I,J)*HUxyHVyy+V(I,J)*HUxyHVyy_y) &
                +ken5 & 
                   *(U(I,J)*UxxVxy+V(I,J)*UxyVyy) &
                +ken4  &
                   *(Uy(I,J)*UxxVxy+U(I,J)*UxxVxy_y+Vy(I,J)*UxyVyy+V(I,J)*UxyVyy_y) &
                +Beta_2*Beta_2*(DUx(I,J)+DVy(I,J)+reta*Beta_2*(Ux(I,J)+Vy(I,J)))  &
                   *(HUxyHVyy+ETAy(I,J)*Beta_2*(Ux(I,J)+Vy(I,J))+reta*Beta_2*UxyVyy)

        omega_0=Vx(I,J)-Uy(I,J)

! ykchoi (15. 08. 06.)
! omega_1 term is modified as Shi et al. (2012) publication
        omega_1=( b2*rhx + Beta_1*ETAx(I,J) )*Beta_2*( HUxyHVyy + (b2*rh+Beta_1*reta)*Beta_2*UxyVyy )  &
              - ( b2*rhy + Beta_1*ETAy(I,J) )*Beta_2*( HUxxHVxy + (b2*rh+Beta_1*reta)*Beta_2*UxxVxy )

	omega(I,J)=omega_0+omega_1

        IF(OUT_VORmax) THEN
        IF(abs(omega(I,J)).GT.VorticityMax(I,J)) THEN
        VorticityMax(I,J)=omega(I,J)
        ENDIF
        ENDIF

       ken1=((Beta_1-1.0_SP/2.0_SP)*(reta+rh)*Beta_2)
       ken2=(1.0_SP/3.0_SP-Beta_1+0.5_SP*Beta_1*Beta_1)*rh*rh*Beta_2*Beta_2  &
               + (1.0_SP/6.0_SP-Beta_1+Beta_1*Beta_1)*rh*reta*Beta_2 &
               +(1.0_SP/2.0_SP*Beta_1*Beta_1-1.0_SP/6.0_SP)*reta*reta*Beta_2*Beta_2

       U3(I,J)=-V(I,J)*omega_1 - omega_0 &
                 *(ken1*HUxyHVyy &
                   +ken2*UxyVyy)
! ykchoi

	 V3(I,J) = U(I,J)*omega_1 + omega_0 &
	            *(ken1*HUxxHVxy &
                    +ken2*UxxVxy)
! ykchoi
       ENDDO
       ENDDO            

     ENDIF  


END SUBROUTINE CAL_DISPERSION

