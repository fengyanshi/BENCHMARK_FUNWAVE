        PROGRAM MAKE_INITIAL
        USE GLOBAL

        USE INPUT_UTIL
        IMPLICIT NONE
        INTEGER :: I,J,K,line
        CHARACTER(LEN=80) File_Name
        INTEGER :: m,n,m1,n1,l,l1
        REAL,DIMENSION(:,:),ALLOCATABLE :: var2D
        REAL,DIMENSION(:,:,:),ALLOCATABLE :: var3D
        REAL :: xx

     ! read from input.txt
       FILE_NAME='input.txt'

     ! dimension                                             
       CALL GET_INTEGER_VAL(Mglob,FILE_NAME,'Mglob',line)
       CALL GET_INTEGER_VAL(Nglob,FILE_NAME,'Nglob',line)
       CALL GET_INTEGER_VAL(Kglob,FILE_NAME,'Kglob',line)   

       m=Mglob
       n=Nglob
       l=Kglob
       ALLOCATE(var2D(m,n),var3D(m,n,l))


      print*,'dimensions',m,n,l
         do j=1,n
         do i=1,m
           xx=(i-1.0)/(m-1.0)*20.0
           var2D(i,j) = 0.1*COS(3.1415926*xx/10.0)
         enddo
         enddo
         
! eta 

100      format(1000f12.3)

         open(2,file='../initial/eta0.txt')

         do j=1,n
          write(2,100)(var2D(i,j),i=1,m)
         enddo
         close(2)

         do k=1,l
         do j=1,n
         do i=1,m
           var3D(i,j,k)=0.0
         enddo
         enddo
         enddo

!     

         open(2,file='../initial/u0.txt')

!         do k=1,l
         k=5
         do j=1,n
          write(2,100)(var3D(i,j,k),i=1,m)
         enddo
!         enddo
         close(2)
         open(2,file='../initial/v0.txt')

!         do k=1,l
         k=5
         do j=1,n
          write(2,100)(var3D(i,j,k),i=1,m)
         enddo
!         enddo
         close(2)
 

         END

