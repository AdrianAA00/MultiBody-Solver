program main
  implicit none
  ! --------------------------------------------------------------------- !
  ! Variables declaration
  real*8, parameter :: pi=acos(-1d0)

  !Initial conditions
  real*8, parameter :: psi01 = 0d0* (pi/180d0), theta01 = 0d0* (pi/180d0), phi01 = 0d0* (pi/180d0)     !Master
  real*8, parameter :: dpsi01 = 0d0* (pi/180d0), dtheta01 = 0d0* (pi/180d0), dphi01 = 0d0* (pi/180d0)
  real*8, parameter :: psi02 = 0d0* (pi/180d0), theta02 = 90d0* (pi/180d0), phi02 = 0d0* (pi/180d0)     !Slave
  real*8, parameter :: dpsi02 = 0d0* (pi/180d0), dtheta02 = 0d0* (pi/180d0), dphi02 = 0d0* (pi/180d0)

  !Numerical solution
  integer :: iter_max=1000, iter
  real*8 :: delta_t=0.001, gamma=0.5, beta=0.25, tau=0.99 !Relaxation factor. First order: gamma = 0.75 and beta = 0.390625

  real*8 :: M(22,22), b(22), aux(22,22)
  real*8 :: x(22), x_new(22)  !Solution to the system
  integer :: pivot(21)

  !Physical problem parameters
  real*8, parameter :: g=9

  real*8, parameter :: l1=1d0, I1=1d0, mass1 = 1d0 !Master
  real*8, parameter :: l2=1d0, I2=1d0, mass2 = 1d0 !Slave

  !Inertia matrices
  real*8 :: mRR1(3,3)=reshape(    (/ 1d0, 0d0, 0d0, &
                                    0d0, 1d0, 0d0, &
                                    0d0, 0d0, 1d0 /), &
                                (/3,3/), order = (/2,1/) )

  real*8 :: m001(3,3)=reshape(    (/ 1d0, 0d0, 0d0, &
                                    0d0, 1d0, 0d0, &
                                    0d0, 0d0, 1d0 /), &
                                (/3,3/), order = (/2,1/) )

  real*8 :: mRR2(3,3)=reshape(    (/ 1d0, 0d0, 0d0, &
                                    0d0, 1d0, 0d0, &
                                    0d0, 0d0, 1d0 /), &
                                (/3,3/), order = (/2,1/) )

  real*8 :: m002(3,3)=reshape(    (/ 1d0, 0d0, 0d0, &
                                    0d0, 1d0, 0d0, &
                                    0d0, 0d0, 1d0 /), &
                                (/3,3/), order = (/2,1/) )

  !Forces and moments
  real*8 :: Q_e1(7)=reshape(  (/ 0d0, 0d0, -g*mass1, 0d0, 0d0, 0d0, 0d0/), (/7/) )  !Master
  real*8 :: Q_e2(7)=reshape(  (/ 0d0, 0d0, -g*mass2, 0d0, 0d0, 0d0, 0d0/), (/7/) )  !Slave

  !Constrains
  real*8, parameter :: x_p1(3)=reshape(  (/ 0d0 , 0d0, l1 /), (/3/) )        !Master
  real*8, parameter :: ddr_p1(3)=reshape(  (/ 0d0 , 0d0, 0d0 /), (/3/) ) 
  real*8, parameter :: x_p2(3)=reshape(  (/ 0d0 , 0d0, l2 /), (/3/) )        !Slave
  real*8, parameter :: ddr_p2(3)=reshape(  (/ 0d0 , 0d0, 0d0 /), (/3/) ) 

  !Counts
  integer :: print
  integer :: iter2
  integer :: Ierr

  !Initial conditions vectorization
  real*8 :: EuAng1(3)=reshape(  (/ psi01,  theta01, phi01 /), (/3/) ) !Master
  real*8 :: EuVel1(3)=reshape(  (/ dpsi01, dtheta01 , dphi01  /), (/3/) )
  real*8 :: EuAng2(3)=reshape(  (/  psi02 , theta02, phi02   /), (/3/) ) !Slave
  real*8 :: EuVel2(3)=reshape(  (/ dpsi02, dtheta02, dphi02   /), (/3/) )
  
  real*8 :: r1(3), dr1(3) !Inertial frame position and coordinate, mass centre
  real*8 :: e1(4), de1(4) !Euler parameters and derivatives
  real*8 :: r2(3), dr2(3) !Inertial frame position and coordinate, mass centre
  real*8 :: e2(4), de2(4) !Euler parameters and derivatives

  real*8 :: q(22), dq(22), q_temp(22), dq_temp(22), q_temp_new(22), dq_temp_new(22) !Vector containing all coordinates
  real*8 :: q_res, dq_res  !Residuals

  real*8 :: ddc1, ddc2  !Euler velocity constrain

  real*8, parameter :: TOLERANCE = 10.0*epsilon(1.0_8)  !Implicit method tolerance

  ! --------------------------------------------------------------------- !
  ! Initialize initial condition vectors

  !Master
  e1 = EulerParam(EuAng1)
  de1 = dEulerParam(EuAng1,EuVel1)
  r1 = MATMUL(fA(e1(1),e1(2),e1(3),e1(4)),-x_P1)
  dr1 =0d0

  !Slave
  e2 = EulerParam(EuAng2)
  de2 = dEulerParam(EuAng2,EuVel2)
  r2 = MATMUL(fA(e2(1),e2(2),e2(3),e2(4)),-x_P2) + r1
  dr2 =0d0

  !All variables
  q = 0d0
  dq = 0d0
  q(1:3) = r1   !Master
  q(4:7) = e1
  dq(1:3) = dr1
  dq(4:7) = de1 
  q(12:14) = r2   !Slave
  q(15:18) = e2
  dq(12:14) = dr2
  dq(15:18) = de2

  !Euler parameters velocity correction
  ddc1 = 0d0
  ddc2 = 0d0

  !Inertia and mass properties scaling
  mRR1 = mRR1*mass1   !Master
  m001 = m001*I1
  mRR2 = mRR2*mass2   !Slave
  m002 = m001*I2

  ! --------------------------------------------------------------------- !
  ! Verification of matrices and solver

  write(*,*) "________________Master________________"
  write(*,*) "Euler Initial Angles =", EuAng1
  write(*,*) "Euler Initial Parameter =", e1
  write(*,*) "Euler Initial Angles velocities=", EuVel1
  write(*,*) "Euler Initial Parameter velocities =", de1

  write(*,*) "________________Slave________________"
  write(*,*) "Euler Initial Angles =", EuAng2
  write(*,*) "Euler Initial Parameter =", e2
  write(*,*) "Euler Initial Angles velocities=", EuVel2
  write(*,*) "Euler Initial Parameter velocities =", de2

  M(1:11,1:11) = fM(x_p1,e1,mRR1,m001)
  M(12:22,12:22) = fM(x_p2,e2,mRR2,m002)
  M(1,20) = -1d0
  M(2,21) = -1d0
  M(3,22) = -1d0
  M(20,1) = -1d0
  M(21,2) = -1d0
  M(22,3) = -1d0

  aux = M
  write(*,*) "________M_________"
  write(*,'(22F8.2)') Transpose(aux(:,:))
  !write(*,*) "Det =", DET(aux)
  
  b(1:11) = fb(x_p1,e1,de1,Q_e1,ddc1,m001,ddr_p1)
  b(12:22) = fb(x_p2,e2,de2,Q_e2,ddc2,m002,ddr_p2)
  write(*,*) "________b_________"
  write(*,'(22F8.2)') b(:)
  

  call LU_Fact(M,pivot,Ierr)
  write(*,*) "________LU_________"
  write(*,'(22F8.2)') Transpose(M(:,:))

  b = pivoting(b,pivot)
  call solve(M, b, x)

  
  write(*,*) "______verify solver______" 
  write(*,'(22F8.5)') MATMUL(aux,x) - b

  write(*,*) "________solution________" 
  write(*,'(22F9.5)') x

  ! --------------------------------------------------------------------- !
  ! Numerical solution

  !Opening file to write results
  open (action='write', file="data_numerical.txt", newunit=print, status='replace')

  x_new = 0

  do iter = 1, iter_max
    q_temp = q
    dq_temp = dq
    x = x_new
    
    iter2 = 0
    q_res = 1.0; dq_res = 1.0
    do while( ( dq_res > TOLERANCE .or. q_res > TOLERANCE ) .and. iter2 < 100)
      iter2 = iter2+1
      ! Normalize
      !Solving
      M = 0d0
      M(1:11,1:11) = fM(x_p1,q_temp(4:7),mRR1,m001)
      M(12:22,12:22) = fM(x_p2,q_temp(15:18),mRR2,m002)
      M(1,20) = -1d0
      M(2,21) = -1d0
      M(3,22) = -1d0
      M(20,1) = -1d0
      M(21,2) = -1d0
      M(22,3) = -1d0
      aux=M
      !write(*,*) "Det =", DET(aux)
      
      ddc1 = ddc_new(q(4:7),dq(4:7),((1-gamma)*x(4:7) + gamma*x_new(4:7)),delta_t,gamma,beta)
      ddc2 = ddc_new(q(15:18),dq(15:18),((1-gamma)*x(15:18) + gamma*x_new(15:18)),delta_t,gamma,beta) 
      b(1:11) = fb(x_p1,q_temp(4:7),dq_temp(4:7),Q_e1,ddc1,m001,ddr_p1)
      b(12:22) = fb(x_p2,q_temp(15:18),dq_temp(15:18),Q_e2,ddc2,m002,ddr_p2)

      call LU_Fact(M,pivot,Ierr)
      if (Ierr /= 0) then
        write(*,*) "Error in LU"
        stop
      end if

      b = pivoting(b,pivot)
      call solve(M, b, x_new)

      !Updating values
      q_temp_new = q + delta_t*(dq) + (delta_t**2)*((1-2*beta)*x + 2*beta*x_new)/2
      dq_temp_new = dq + delta_t*((1-gamma)*x + gamma*x_new)
      q_temp = (1d0-tau)*q_temp + tau*q_temp_new 
      dq_temp = (1d0-tau)*dq_temp + tau*dq_temp_new 
      q_temp(4:7) = q_temp(4:7) / sqrt(sum(q_temp(4:7)**2))

      !write (*, *)  "iter", iter, "=", x_new
      ! write(*,*) "   ",iter2, "residual = ", q_res , dq_res
    end do

    q = q_temp
    dq = dq_temp
    EuAng1 = EulerAng(q(4:7))
    EuAng2 = EulerAng(q(15:18))
    !write(*,*) "theta*theta =", DOT_PRODUCT(q_temp(4:7),q_temp(4:7))
    !write(*,*) "theta*theta_dot =", 2d0*DOT_PRODUCT(q_temp(4:7),dq_temp(4:7))
    !write (*, *)  "q",  "=", q
    !write (*, *)  "    ",  q(1:3)
    !write (*,'(11E10.2)') q
    write (print, *) q(13) ,q(14) !pi/2d0 + atan2(q(3),q(2)), pi/2d0 + atan2((q(14)-q(3)),(q(13)-q(2))) !atan(q(3)/q(2)), atan((q(14)-q(3))/(q(13)-q(2))) !EuAng1(2), EuAng2(2)

  end do

  close (print)

  contains

  subroutine analytical (theta0,N,t,g,l)
    implicit none
      integer, intent(in) :: N
      real*8, intent (in) :: theta0, t, g, l
      real*8, allocatable :: theta(:)
      integer :: i,fu

      allocate (theta(N))

      do i=1,N
        theta(i)=abs(theta0*cos(i*t*sqrt(g/l)))
      end do

      open (action='write', file="data_analytical.txt", newunit=fu, status='replace')

      do i = 1, N
          write (fu, *) i, theta(i)
      end do

      close (fu)

  end subroutine analytical 

  function fG_ (e0,e1,e2,e3)
    implicit none
      real*8, intent (in) :: e0, e1, e2, e3
      real*8 :: fG_(3,4)

      fG_=reshape(    (/ -e1, e0, e3, -e2,&
                         -e2, -e3, e0, e1,&
                         -e3, e2, -e1, e0   /), &
                           (/3,4/), order = (/ 2, 1 /) )

      fG_ = 2*fG_

  end function fG_ 

  function fG (e0,e1,e2,e3)
    implicit none
      real*8, intent (in) :: e0, e1, e2, e3
      real*8 :: fG(3,4)

      fG=reshape(    (/ -e1, e0, -e3, e2,&
                        -e2, e3, e0, -e1,&
                        -e3, -e2, e1, e0  /), &
                           (/3,4/), order = (/2,1/) )

      fG = 2*fG

  end function fG 

  function dfG_ (de0,de1,de2,de3)
    implicit none
      real*8, intent (in) :: de0, de1, de2, de3
      real*8 :: dfG_(3,4)

      dfG_=reshape(    (/ -de1, de0, de3, -de2,&
                          -de2, -de3, de0,  de1,&
                          -de3, de2, -de1, de0/), &
                           (/3,4/), order = (/ 2, 1 /) )

      dfG_ = 2*dfG_

  end function dfG_

  function fA (e0,e1,e2,e3)
    implicit none
      real*8, intent (in) :: e0, e1, e2, e3
      real*8 :: fA(3,3)

      fA(1,1) = 2d0*(e0**2 + e1**2) - 1d0
      fA(2,1) = 2d0*(e1*e2 + e0*e3)
      fA(3,1) = 2d0*(e1*e3 - e0*e2)
      fA(2,2) = 2d0*(e0**2 + e2**2) - 1d0
      fA(3,2) = 2d0*(e2*e3 + e0*e1)
      fA(3,3) = 2d0*(e0**2 + e3**2) - 1d0
      fA(1,2) = 2d0*(e1*e2 - e0*e3)
      fA(1,3) = 2d0*(e1*e3 + e0*e2)
      fA(2,3) = 2d0*(e2*e3 - e0*e1)

  end function fA 

  function dfA (e0,e1,e2,e3,de0,de1,de2,de3)
    implicit none
      real*8, intent (in) :: e0, e1, e2, e3,de0,de1,de2,de3
      real*8 :: dfA(3,3)

      dfA(1,1) = 4d0*(e0*de0 + e1*de1)
      dfA(2,1) = 2d0*(de1*e2 + e1*de2 + de0*e3 + e0*de3)
      dfA(3,1) = 2d0*(de1*e3 + e1*de3 - e0*de2 - de0*e2)
      dfA(2,2) = 4d0*(e0*de0 + e2*de2)
      dfA(3,2) = 2d0*(de2*e3 + e2*de3 + de0*e1 + e0*de1)
      dfA(3,3) = 4d0*(e0*de0 + e3*de3)
      dfA(1,2) = 2d0*(de1*e2 + e1*de2 - e0*de3 - de0*e3)
      dfA(1,3) = 2d0*(de1*e3 + e1*de3 + de0*e2 + e0*de2)
      dfA(2,3) = 2d0*(de2*e3 + e2*de3 - de0*e1 - e0*de1)

  end function dfA 

  function fu (x1,x2,x3)
    implicit none
      real*8, intent (in) :: x1, x2 ,x3
      real*8 :: fu(3,3)

      fu(1,1) = 0d0
      fu(2,1) = x3
      fu(3,1) = -x2
      fu(2,2) = 0d0
      fu(3,2) = x1
      fu(3,3) = 0d0
      fu(1,2) = -x3
      fu(1,3) = x2
      fu(2,3) = -x1

  end function fu 
  
  function fM(x,e, mRR, m00)
    implicit none
      real*8, intent (in) :: x(3), e(4), mRR(3,3), m00(3,3)
      real*8 :: e0, e1, e2, e3
      real*8 :: fM(11,11), u(3,3), A(3,3), G(3,4), G_(3,4)

      e0 = e(1)
      e1 = e(2)
      e2 = e(3)
      e3 = e(4)

      fM = 0d0
      G_ = fG_ (e0,e1,e2,e3)
      A = fA (e0,e1,e2,e3)
      u = fu (x(1),x(2),x(3))
      
      fM(1:3,1:3) = mRR
      fM(4:7,4:7) = MATMUL(MATMUL(TRANSPOSE(G_),m00),G_)
      fM(8,4) = e0
      fM(8,5) = e1
      fM(8,6) = e2
      fM(8,7) = e3
      fM(4,8) = e0
      fM(5,8) = e1
      fM(6,8) = e2
      fM(7,8) = e3
      fM(9,1) = 1d0
      fM(10,2) = 1d0
      fM(11,3) = 1d0
      fM(1,9) = 1d0
      fM(2,10) = 1d0
      fM(3,11) = 1d0
      fM(9:11,4:7) = -MATMUL(MATMUL(A,u),G_)
      fM(4:7,9:11) = -TRANSPOSE(MATMUL(MATMUL(A,u),G_))

  end function fM

  function fb(x,e,de,Q,ddc,m00,ddr_p)
    implicit none
      real*8, intent (in) :: x(3), e(4), de(4)
      real*8 :: e0, e1, e2, e3, de0, de1, de2, de3, ddc
      real*8 :: fb(11), u(3,3), dA(3,3), G_(3,4), Q(7), dG_(3,4), m00(3,3), ddr_p(3)

      e0 = e(1)
      e1 = e(2)
      e2 = e(3)
      e3 = e(4)
      de0 = de(1)
      de1 = de(2)
      de2 = de(3)
      de3 = de(4)

      fb = 0
      G_ = fG_ (e0,e1,e2,e3)
      dG_ = dfG_ (de0,de1,de2,de3)
      dA = dfA (e0,e1,e2,e3,de0,de1,de2,de3)
      u = fu (x(1),x(2),x(3))

      fb(1:7) = Q
      fb(4:7) = fb(4:7) - 2d0*MATMUL(MATMUL(MATMUL(TRANSPOSE(dG_),m00),G_), de)
      fb(8) = -DOT_PRODUCT(de,de) + 0.5*ddc
      fb(9:11) = MATMUL(MATMUL(MATMUL(dA,u),G_), de) + ddr_p 
!
  end function fb

  function EulerParam(e_a)
    implicit none
      real*8, intent (in) :: e_a(3)
      real*8 :: EulerParam(4)
      
      EulerParam(1) = cos(e_a(2)/2d0)*cos((e_a(1) + e_a(3))/2d0)
      EulerParam(2) = sin(e_a(2)/2d0)*cos((e_a(1) - e_a(3))/2d0)
      EulerParam(3) = sin(e_a(2)/2d0)*sin((e_a(1) - e_a(3))/2d0)
      EulerParam(4) = cos(e_a(2)/2d0)*sin((e_a(1) + e_a(3))/2d0)

  end function EulerParam

  function dEulerParam(e_a, de_a)
    implicit none
      real*8, intent (in) :: e_a(3), de_a(3)
      real*8 :: dEulerParam(4)
      
      dEulerParam(1) = -sin(e_a(2)/2)*cos((e_a(1) + e_a(3))/2)*de_a(2)/2 - &
      cos(e_a(2)/2)*sin((e_a(1) + e_a(3))/2)*(de_a(1) + de_a(3))/2
      dEulerParam(2) = cos(e_a(2)/2)*cos((e_a(1) - e_a(3))/2)*de_a(2)/2 - &
      sin(e_a(2)/2)*sin((e_a(1) - e_a(3))/2)*(de_a(1) - de_a(3))/2
      dEulerParam(3) = cos(e_a(2)/2)*sin((e_a(1) - e_a(3))/2)*de_a(2)/2 + &
      sin(e_a(2)/2)*cos((e_a(1) - e_a(3))/2)*(de_a(1) - de_a(3))/2
      dEulerParam(4) = -sin(e_a(2)/2)*sin((e_a(1) + e_a(3))/2)*de_a(2)/2 + &
      cos(e_a(2)/2)*cos((e_a(1) + e_a(3))/2)*(de_a(1) + de_a(3))/2

  end function dEulerParam

  function EulerAng(e)
    implicit none
      real*8, intent (in) :: e(4)
      real*8 :: EulerAng(3), A(3,3)
      
      A = fA (e(1),e(2),e(3),e(4))
      EulerAng(1) = datan(-A(1,3)/A(2,3))
      EulerAng(2) = datan(sqrt(abs(1d0-A(3,3)**2d0))/ A(3,3))
      EulerAng(3) = datan(A(3,1)/ A(3,2))

  end function EulerAng

  subroutine LU_Fact(A,Pivot,Ierr)
    IMPLICIT NONE
    !
    ! The arguments:
    !
    ! The input matrix, at exit
    ! will hold the LU factorization
    REAL*8, INTENT(INOUT) :: A(:,:)
    ! Vector of permutations
    INTEGER, INTENT(OUT)   :: Pivot(:)
    ! Singularity indicator, = 0 if A nonsingular,
    ! and = column number j if the first zero
    ! pivot was detected in column j
    INTEGER, INTENT(OUT)   :: Ierr
    !
    ! The local variables:
    !
    LOGICAL :: singular        ! Singularity flag 
    INTEGER :: i,j,k,n         ! DO loop variables
    INTEGER, ALLOCATABLE :: p(:) ! pivot location
    REAL*8, ALLOCATABLE :: Tmp(:)  ! Temporary row
    REAL*8  :: uround            ! rounding unit
    !
    ! Check if the argument is a square matrix
    IF( size(A,1).NE.size(A,2) ) THEN
      PRINT*,"Error in Factorize: A must be square"
      RETURN
    END IF
    n=SIZE(A,1) ! the dimension of the matrix

    ALLOCATE(Tmp(n),stat = k)
    IF (k.NE.0) THEN
      PRINT*,"Error in Factorize: cannot allocate Tmp"; RETURN
    END IF
    !
    ALLOCATE(p(n),stat = k)
    IF (k.NE.0) THEN
      PRINT*,"Error in Factorize: cannot allocate p"; RETURN
    END IF
    !
    Ierr = 0                      ! reset error indicator
    singular = .FALSE.            ! reset singularity flag
    uround = 1.0E-7               ! unit roundoff, set it to
                                  ! 1.0D-14 for double precision
    !
    DO j=1,n-1  ! columns 1:n-1
    !
      p=MAXLOC(ABS(A(j:n,j)))+j-1 ! Look for pivot, A(p(1),j)
      IF (p(1).NE.j) THEN         ! If pivot is non-diagonal
         Tmp(:)    = A(j,:)       ! permute rows j and p(1)
         A(j,:)    = A(p(1),:)
         A(p(1),:) = Tmp(:)
         Pivot(j) = p(1)          ! Save the pivot position
      ELSE
         Pivot(j) = j             ! Save the pivot position
      END IF
    !
    ! If A(j,j)=0 then the matrix is singular
    ! uround is the rounding unit (machine epsilon)
      IF ( ABS(A(j,j)) .LT. uround  ) THEN
        Ierr = j                  ! Singularity Position
        singular = .TRUE.         ! Singularity Flag
        exit                      ! the `DO j' loop
      END IF
    !
      DO i=j+1,n                  ! rows to be processed
         A(i,j) = A(i,j)/A(j,j)   ! Store the multiplier
         A(i,j+1:n) = A(i,j+1:n) - A(i,j)*A(j,j+1:n)
      END DO          ! Row loop
    !
    END DO            ! Column loop
    !
    ! If A(n,n)=0 then the matrix is singular
    ! uround is the rounding unit (machine epsilon)
    IF ( abs(A(n,n)) .LT. uround  ) THEN
        Ierr = n                  ! Singularity Flag
        singular = .TRUE.
    END IF
    !  
    IF (allocated(Tmp)) DEALLOCATE(Tmp)
    !
    IF (.NOT. singular) THEN
       Ierr = 0
    END IF

  end subroutine

  function pivoting(b,pivot)
    implicit none
      real*8, intent (in) :: b(:)
      integer, intent (in) :: pivot(:)
      real*8, allocatable :: pivoting(:)
      real*8, allocatable :: aux(:)
      integer :: n,i

      n=SIZE(b)
      allocate(aux(n-1))
      allocate(pivoting(n))

      pivoting = b
      DO i=1,n-1 
         aux = pivoting                 
         pivoting(i) = b(pivot(i))
         pivoting(pivot(i)) = aux(i)
      END DO  

  end function pivoting

  subroutine solve(a, b, x)
    implicit none
    real*8, intent(inout)  :: a(:, :)
    real*8, intent(in)  :: b(:)
    real*8, intent(out) :: x(:)
    integer :: n, i, j
    real*8 :: s
    real*8, allocatable :: y(:), L(:,:), U(:,:)

    n = size(x)

    allocate(y(n),L(n,n),U(n,n))

    do i = 1, n
      s = sum((/(a(i, j) * y(j), j=1,i-1)/))
      y(i) = b(i) - s
    end do

    do i = n, 1, -1
      s = sum((/(a(i, j) * x(j), j=i+1,n)/))
      x(i) = (y(i) - s) / a(i, i)
    end do

    deallocate(y)

  end subroutine solve

  real*8 function DET(aa)
    real*8 aa(:,:)
    real*8 tmp,c(size(aa,dim=1),size(aa,dim=2))
    real*8 max
    integer i,j,k,l,m,n,num(size(aa,dim=1))
      n=size(aa,dim=1)
      det=1.	
      do k=1,n
        max=aa(k,k);num(k)=k;
        do i=k+1,n 
          if(abs(max)<abs(aa(i,k))) then
            max=aa(i,k)
            num(k)=i
          endif
        enddo
        if (num(k)/=k) then
          do l=k,n 
            tmp=aa(k,l)
            aa(k,l)=aa(num(k),l)
            aa(num(k),l)=tmp
          enddo
          det=-1.*det
        endif
        do m=k+1,n
          c(m,k)=aa(m,k)/aa(k,k)
          do l=k,n 
            aa(m,l)=aa(m,l)-c(m,k)*aa(k,l)
          enddo
        enddo !There we made matrix triangular!	
      enddo

      do i=1,n
      det=det*aa(i,i)
      enddo
      return

  end function

  function ddc_new(e,de,dde,t,gamma,beta)
    implicit none
      real*8, intent (in) :: e(:), de(:), dde(:)
      real*8 :: ddc_new, c0, c1, c2, d0, d1, cn, cn1, cn2, t, gamma, beta, A, B, C

      c0 = t
      c1 = 0.5*t**2d0*(1d0-2d0*beta)
      c2 = t**2d0*beta
      d0 = t*(1d0-gamma)
      d1 = t*gamma

      cn = DOT_PRODUCT(e,e)
      cn1 = 2d0*DOT_PRODUCT(de,e)
      cn2 = 2d0*(DOT_PRODUCT(de,de) + DOT_PRODUCT(dde,e))
      !write(*,*) cn, cn1, cn2

      A = 1d0 - cn - 2d0*c0*cn1 - (c1 + c0*d0)*cn2
      B = c2 + c0*d1 + c1 - (d0+d1)*c2/d1
      C = (c2/d1)*(cn1 + d0*cn2)
      !write(*,*) A,B,C

      !ddc_new = (A+C)/B
      ddc_new = -(cn1/(gamma*t)+(1d0-gamma)*cn2/gamma)

      !write(*,*) ddc_new

  end function ddc_new

end
