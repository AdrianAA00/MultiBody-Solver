program main
  implicit none
  
  integer ::  iter_max
  real*8    :: delta_t

  integer  :: NumberBodies, size_M, Ierr, count, count2

  integer, allocatable  :: constrains(:,:), constrain_coupling(:,:)

  real*8 , allocatable  :: initial_conditions(:,:)

  real*8 , allocatable :: fixed_coordinates(:,:)

  real*8 , allocatable :: inertia_tensors(:,:)

  real*8 , allocatable :: mass(:,:)

  real*8 , allocatable :: Q_g(:,:)
  
  real*8 , allocatable :: M(:,:), b(:), q(:,:), dq(:,:)

  real*8 , allocatable :: x_p(:,:,:)

  integer, allocatable :: pivot(:)
  integer :: iter, iter2, print
  real*8 :: q_res, dq_res

  real*8, allocatable :: x_new(:), auxM(:,:)
  real*8, allocatable :: q_temp(:,:), q_temp_new(:,:), dq_temp(:,:), dq_temp_new(:,:), ddq(:,:), ddq_new(:,:)

  real*8, allocatable :: ddconstrains(:,:), ddconstrains_new(:,:)
  real*8, allocatable :: q_constrains(:,:), dq_constrains(:,:), q_constrains_temp(:,:), dq_constrains_temp(:,:)

  real*8, parameter :: beta = 0.25, gamma = 0.5, tau = 0.99   !Newmark parameters 
  real*8, parameter :: TOLERANCE = 10.0*epsilon(1.0_8)        !Tolerance implicit scheme

  real*8  :: ddc  !Corrections velocity 

  ! --------------------------------------------------------------------- !
  ! Verification of matrices and solver

  !Reading values from imput file
  call read_initial (iter_max, delta_t, NumberBodies, size_M, constrains, initial_conditions, &
  fixed_coordinates, inertia_tensors, mass)

  allocate(pivot(size_M-1), x_new(size_M), auxM(size_M, size_M))

  allocate(q_temp(8,NumberBodies+1), q_temp_new(8,NumberBodies+1), dq_temp(8,NumberBodies+1))
  allocate(dq_temp_new(8,NumberBodies+1),ddq(8,NumberBodies+1),ddq_new(8,NumberBodies+1))

  allocate(ddconstrains(9,NumberBodies))
  allocate(q_constrains(9,NumberBodies), dq_constrains(9,NumberBodies))
  allocate(q_constrains_temp(9,NumberBodies), dq_constrains_temp(9,NumberBodies))

  !Setting initial conditions
  call initial_condit (initial_conditions, NumberBodies, q, dq)

  !Setting relationships between bodies and reference points
  call x_position_body (NumberBodies, q, fixed_coordinates)  

  !Setting constrains couplings
  call cons_coupling (NumberBodies, fixed_coordinates, constrain_coupling) 
  
  !Defining matrix M
  call fM_general(M, x_p, q, mass, inertia_tensors, constrains, NumberBodies, size_M, constrain_coupling)
  auxM = M
  write(*,*)  "_______________________Matrix M_______________________"
  write(*,"(22F8.3)") Transpose(M)

  !Defining vector b
  call read_iter(Q_g, ddconstrains_new, NumberBodies,q)
  call fb_general (b,x_p,q,dq,NumberBodies, size_M, Q_g, ddconstrains_new, constrain_coupling)
  write(*,*)  "_______________________Vector b_______________________"
  write(*,"(22F8.3)") b
  
  !Solving
  call LU_Fact(M,pivot,Ierr)

  b = pivoting(b,pivot)
  call solve(M, b, x_new)
  call f_ddq(x_new,ddq, NumberBodies, constrains, initial_conditions)

  write(*,*) "______verify solver______" 
  write(*,'(22F8.5)') MATMUL(auxM,x_new) - b

  write(*,*) "________solution________" 
  write(*,'(22F9.5)') x_new

  ! --------------------------------------------------------------------- !
  ! Numerical solution

  !Opening file to write results
  open (action='write', file="data_numerical_general.txt", newunit=print, status='replace')

  ddq_new = 0

  do iter = 1, iter_max
    !Coordinates update
    q_temp = q
    dq_temp = dq
    ddq = ddq_new
    q_constrains_temp = q_constrains
    dq_constrains_temp = q_constrains
    ddconstrains = ddconstrains_new

    iter2 = 0
    q_res = 1d0; dq_res = 1d0
    do while( ( dq_res > TOLERANCE .or. q_res > TOLERANCE ) .and. iter2 < 1)
      iter2 = iter2+1

      !Defining matrix M
      M = 0d0
      call fM_general(M, x_p, q_temp, mass, inertia_tensors, constrains, NumberBodies, size_M, constrain_coupling)
      
      !Defining matrix b
      b = 0d0
      call read_iter(Q_g, ddconstrains_new, NumberBodies,q)
      call fb_general (b,x_p,q_temp,dq_temp,NumberBodies, size_M, Q_g, ddconstrains_new, constrain_coupling)

      !Solving
      call LU_Fact(M,pivot,Ierr)
      if (Ierr /= 0) then
        write(*,*) "Error in LU"
        stop
      end if

      b = pivoting(b,pivot)

      call solve(M, b, x_new)
      call f_ddq(x_new,ddq_new, NumberBodies, constrains, initial_conditions)

      !Updating values, Newmark
      q_temp_new(2:,:) = q(2:,:) + delta_t*(dq(2:,:)) + (delta_t**2)*((1-2*beta)*ddq(2:,:)+2*beta*ddq_new(2:,:))/2
      dq_temp_new(2:,:) = dq(2:,:) + delta_t*((1-gamma)*ddq(2:,:) + gamma*ddq_new(2:,:))
      q_temp(2:,:) = (1d0-tau)*q_temp(2:,:) + tau*q_temp_new(2:,:) 
      dq_temp(2:,:) = (1d0-tau)*dq_temp(2:,:) + tau*dq_temp_new(2:,:) 

      !Correction for norm of euler parameters
      do count = 2, NumberBodies + 1
        q_temp(5:8,count) = q_temp(5:8,count) / sqrt(sum(q_temp(5:8,count)**2))
      end do

      !Constrain correction, explicit Newmark integration
      q_constrains_temp(3:,:) = q_constrains(3:,:) + delta_t*(dq_constrains(3:,:)) + &
      (delta_t**2)*((1-2*beta)*ddconstrains(3:,:)+2*beta*ddconstrains_new(3:,:))/2
      dq_constrains_temp(3:,:) = dq_constrains(3:,:) + delta_t* &
      ((1-gamma)*ddconstrains(3:,:) + gamma*ddconstrains_new(2:,:))
      call position_correction (q_temp,constrains,x_p,NumberBodies,q_constrains)

      !Residuals
      q_res = sqrt(sum((q_temp_new(2:,:) - q_temp(2:,:))**2))
      dq_res = sqrt(sum((dq_temp_new(2:,:) - dq_temp(2:,:))**2))

      !write(*,*) "   ",iter2, "residual = ", q_res , dq_res
    end do

    q(2:,:) = q_temp(2:,:)
    dq(2:,:) = dq_temp(2:,:)
    q_constrains(3:,:) = q_constrains_temp(3:,:)
    dq_constrains(3:,:) = q_constrains_temp(3:,:)

    !EuAng1 = EulerAng(q(4:7))
    !EuAng2 = EulerAng(q(15:18))

    write (print, *) q(3:4,3) 

  end do

  close (print)

  contains

  subroutine position_correction (q_temp,indices,x_p,NumberBodies,q_constrains)
    implicit none
      real*8, allocatable, intent (inout) :: q_temp(:,:)
      real*8, allocatable, intent (in) :: q_constrains(:,:), x_p(:,:,:)
      integer, intent (in) :: indices(:,:), NumberBodies
      real*8 ::  A(3,3)
      real*8 :: auxq(3) 
      integer :: index, count, count2

      auxq = 0d0

      do count2 = 1, NumberBodies +1  !Body attached to fixed point index
        do count = 1, NumberBodies +1  !Fixed Point index
          if ((constrain_coupling(count,count2) == 1 .AND. count == count2)) then
            A = fA (q_temp(5,count),q_temp(6,count),q_temp(7,count),q_temp(8,count))
            auxq = auxq + q_constrains(3:5,count) - MATMUL(A,x_p(:,count,count))  !Diagonal
          elseif ((constrain_coupling(count,count2) == 1 .AND. count /= count2)) then
            A = fA (q_temp(5,count2),q_temp(6,count2),q_temp(7,count2),q_temp(8,count2))
            auxq = auxq + MATMUL(A,x_p(:,count,count2))  !Diagonal !Off-diagonal
          end if 

          do index = 1,3
            if (indices(2+index,count) == 1) then
              q_temp(1+index,count) = auxq(index)
            end if
          end do

        end do
      end do
      
  end subroutine

  subroutine inertial_general_F (Q_g,F_inert,e)
    implicit none

    real*8, intent(in) :: F_inert(7), e(4)
    real*8, intent(out) :: Q_g(8)
    
    Q_g = 0d0
    Q_g(1:4) = F_inert(1:4)
    Q_g(5) = -(1/2)*(e(2)*F_inert(5)+e(3)*F_inert(6)+e(4)*F_inert(7))
    Q_g(6) = (1/2)*(e(1)*F_inert(5)+e(4)*F_inert(6)-e(3)*F_inert(7))
    Q_g(7) = (1/2)*(-e(4)*F_inert(5)+e(1)*F_inert(6)+e(2)*F_inert(7))
    Q_g(8) = (1/2)*(e(3)*F_inert(5)-e(2)*F_inert(6)+e(1)*F_inert(7))
  end subroutine

  subroutine f_ddq (x_new,ddq, NumberBodies, constrains, initial_conditions)
    implicit none

    real*8, allocatable, intent(in) :: x_new(:), initial_conditions(:,:)
    real*8, allocatable, intent(out) :: ddq(:,:)
    integer, intent(in) :: NumberBodies

    integer :: constrain, constrain_tot, count
    integer, allocatable, intent(in) :: constrains(:,:)

    allocate(ddq(8,NumberBodies+1))
    constrain = 0d0
    constrain_tot = 0d0
    ddq = 0d0
    do count = 2, NumberBodies + 1
      ddq(1,count) = initial_conditions(1,count)
      ddq(2:7,count) = x_new(1+constrain_tot+8*(count-2):7+constrain_tot+8*(count-2))
      constrain = sum(constrains(3:,count-1))
      constrain_tot = sum(constrains(3:,count-1)) + constrain_tot
    end do

  end subroutine


  subroutine fM_general(M, x_p, q, mass, inertia_tensors, constrains, NumberBodies, size_M, constrain_coupling)
    implicit none

    real*8  :: identity(3,3)
    integer :: constrain, constrain_tot, count, count2
    integer , allocatable :: edges_matrix(:), edges_size(:)
    real*8, allocatable, intent(out) :: M(:,:)
    real*8, allocatable, intent(in) :: x_p(:,:,:), q(:,:), mass(:,:), inertia_tensors(:,:)
    integer, allocatable, intent(in) :: constrain_coupling(:,:)
    integer, allocatable, intent(in) :: constrains(:,:)
    integer, intent(in) :: NumberBodies, size_M

    allocate(M(size_M,size_M))

    ! Obtaining matrix M
    M=0d0
    allocate(edges_matrix(NumberBodies))
    allocate(edges_size(NumberBodies))
    identity = 0
    do count =1,3
      identity(count, count) = 1
    end do
    constrain = 0
    constrain_tot = 0
    do count = 1, NumberBodies 
      constrain = sum(constrains(3:,count))
      M((count-1)*8+1+constrain_tot:(count-1)*8+1+constrain_tot+constrain,(count-1)*8+1+constrain_tot:(count-1)*8+1+constrain_tot+ & 
      constrain)=fM(q(5:7,count+1), mass(2,count)*identity, inertia_tensors(2:4,(count-1)*3+1:(count-1)*3+3)) !Mass inertias 
      constrain_tot = sum(constrains(3:,count)) + constrain_tot
      edges_matrix(count) = constrain_tot + count*8
      edges_size(count) = constrain
    end do
    do count2 = 1, NumberBodies +1  !Body attached to fixed point index
      do count = 1, NumberBodies +1  !Fixed Point index
        if ((constrain_coupling(count,count2) == 1) .AND. (count2>1) .AND. (count>1) .AND. (count==count2) ) then
          M(edges_matrix(count-1)-edges_size(count-1)+1:edges_matrix(count-1), edges_matrix(count2-1)-edges_size(count2-1)+1-8: &
          edges_matrix(count2-1)-edges_size(count2-1)-1)  = fMconst(x_p(:,count,count2),q(5:7,count2),constrains(3:,count2-1))            !Constrains diagonal
          M(edges_matrix(count2-1)-edges_size(count2-1)+1-8:edges_matrix(count2-1)-edges_size(count2-1)-1, edges_matrix(count-1)- &
          edges_size(count-1)+1:edges_matrix(count-1))=transpose(fMconst(x_p(:,count,count2),q(5:7,count2),constrains(3:,count2-1)))
        else if ((constrain_coupling(count,count2) == 1) .AND. (count2>1) .AND. (count>1) .AND. (count/=count2)) then
          M(edges_matrix(count-1)-edges_size(count-1)+1:edges_matrix(count-1), edges_matrix(count2-1)-edges_size(count2-1)+1-8: &
          edges_matrix(count2-1)-edges_size(count2-1)-1)  = -fMconst(x_p(:,count,count2),q(5:7,count2),constrains(3:,count2-1))           !Constrains off-diago
          M(edges_matrix(count2-1)-edges_size(count2-1)+1-8:edges_matrix(count2-1)-edges_size(count2-1)-1, edges_matrix(count-1) - &
        edges_size(count-1)+1:edges_matrix(count-1))=-transpose(fMconst(x_p(:,count,count2),q(5:7,count2),constrains(3:,count2-1)))  
        end if 
      end do
    end do

  end subroutine fM_general

  subroutine fb_general (b,x_p,q,dq,NumberBodies,size_M,Q_g,ddconstrains_new, constrain_coupling)
    implicit none
    
    real*8, allocatable, intent(out) :: b(:)
    real*8, allocatable, intent(in) :: x_p(:,:,:), q(:,:), dq(:,:), Q_g(:,:), ddconstrains_new(:,:)
    integer, allocatable, intent(in) :: constrain_coupling(:,:)
    integer , allocatable :: edges_matrix(:), edges_size(:)
    integer, intent(in) :: NumberBodies, size_M
    integer :: count, count2, constrain, constrain_tot

    allocate(b(size_M))
    allocate(edges_matrix(NumberBodies))
    allocate(edges_size(NumberBodies))

    constrain = 0
    constrain_tot = 0
    do count = 1, NumberBodies 
      constrain = sum(constrains(3:,count))
      constrain_tot = sum(constrains(3:,count)) + constrain_tot
      edges_matrix(count) = constrain_tot + count*8
      edges_size(count) = constrain
    end do

    ! Obtaining vector b
    b=0d0
    constrain = 0
    constrain_tot = 0
    do count = 1, NumberBodies 
      constrain = sum(constrains(3:,count))
      b((count-1)*8+1+constrain_tot:(count-1)*8+1+constrain_tot+constrain) = fb(x_p(:,count+1,count+1),q(5:7,count+1) & !Mass inertias + constrains controller
      ,dq(5:7,count+1),Q_g(2:,count), ddc,inertia_tensors(2:4,(count-1)*3+1:(count-1)*3+3))
      constrain_tot = sum(constrains(3:,count)) + constrain_tot
    end do
    do count2 = 1, NumberBodies +1 !Body attached to fixed point index
      do count = 1, NumberBodies +1 !Fixed Point index
        if ((constrain_coupling(count,count2) == 1) .AND. (count2>1) .AND. (count>1) .AND. (count==count2) ) then
          b(edges_matrix(count-1)-edges_size(count-1)+1:edges_matrix(count-1)) &
          = b(edges_matrix(count-1)-edges_size(count-1)+1:edges_matrix(count-1)) + &
          fbconst(x_p(:,count,count2),q(5:7,count2),dq(5:7,count2),constrains(3:,count2-1))           !Constrains diagonal
        else if ((constrain_coupling(count,count2) == 1) .AND. (count2>1) .AND. (count>1) .AND. (count/=count2)) then
          b(edges_matrix(count-1)-edges_size(count-1)+1:edges_matrix(count-1)) &
          = b(edges_matrix(count-1)-edges_size(count-1)+1:edges_matrix(count-1)) - &
          fbconst(x_p(:,count,count2),q(5:7,count2),dq(5:7,count2),constrains(3:,count2-1))           !Constrains off-diago
        end if
      end do
    end do 

    do count = 1, NumberBodies                                                                        !Constrains controller
      b(edges_matrix(count)-edges_size(count)+1:edges_matrix(count)) &
      = b(edges_matrix(count)-edges_size(count)+1:edges_matrix(count)) + &
      fbcontroller(constrains(3:,count), ddconstrains_new(3:,count))
    end do

  end subroutine fb_general


  subroutine cons_coupling (NumberBodies, fixed_coordinates, constrain_coupling)  
    implicit none
    integer :: count
    integer, intent(in) :: NumberBodies
    real*8, allocatable, intent(in) :: fixed_coordinates(:,:)
    integer, allocatable, intent(out) :: constrain_coupling(:,:)

   !Constrain coupling
    allocate(constrain_coupling(NumberBodies+1,NumberBodies+1))
    write(*,*)  "___________________Constrain coupling____________________"
    constrain_coupling = 0
    do count = 1, NumberBodies+1
      constrain_coupling(int(fixed_coordinates(2,count))+1,int(fixed_coordinates(1,count))+1) = 1
      constrain_coupling(int(fixed_coordinates(2,count))+1,int(fixed_coordinates(2,count))+1) = 1
    end do
    write(*,'(3I4)') transpose(constrain_coupling)

  end subroutine cons_coupling

  subroutine x_position (NumberBodies, q, fixed_coordinates)  
    implicit none
    integer :: count, count2
    real*8 :: A(3,3)
    integer, intent(in) :: NumberBodies
    real*8, allocatable, intent(in) :: q(:,:), fixed_coordinates(:,:)

    ! Obtaining x coordinates of fixed points in body frame
    allocate(x_p(3,NumberBodies+1,NumberBodies+1))
    write(*,*)  "_____________Body coordinates fixed points_______________"
    do count = 1, NumberBodies+1   !Reference fixed points
      do count2 = 1, NumberBodies+1    !Centres mass bodies
        A = fA (q(5,count2),q(6,count2),q(7,count2),q(8,count2))
        x_p(:,count,count2) = Matmul(Transpose(A),fixed_coordinates(3:5,count)-q(2:4,count2))
        write(*,'(13F8.3)') x_p(:,count,count2)
      end do
    end do
  end subroutine x_position

  subroutine x_position_body (NumberBodies, q, fixed_coordinates)  
    implicit none
    integer :: count, count2, aux2
    real*8 :: A(3,3)
    real*8, allocatable :: aux(:,:)
    integer, intent(in) :: NumberBodies
    real*8, allocatable, intent(in) :: q(:,:), fixed_coordinates(:,:)

    ! Obtaining x coordinates of fixed points in body frame
    allocate(x_p(3,NumberBodies+1,NumberBodies+1))
    allocate(aux(3,NumberBodies+1))

    ! Converting fixed_coordinates to inertial frame
    write(*,*)  "_____________Inertial coordinates fixed points_______________"
    aux = 0d0
    x_p = 0d0
    do count = 1, NumberBodies + 1
      do count2 = 1, NumberBodies + 1
        aux2 = q(1,count)
        if (int(fixed_coordinates(2, count2)) == aux2) then    !Selecting the body reference frame
          A = fA (q(5,count),q(6,count),q(7,count),q(8,count))
          aux(:, count2) = Matmul(A,fixed_coordinates(3:5,count2)) + q(2:4,count)
          write(*,'(13F8.3)') aux(:,count)
        end if
      end do
    end do

    write(*,*)  "_____________Body coordinates fixed points_______________"
    do count = 1, NumberBodies+1  !Reference points
      do count2 = 1, NumberBodies+1    !Centres mass bodies
        A = fA (q(5,count2),q(6,count2),q(7,count2),q(8,count2))
        x_p(:,count,count2) = Matmul(Transpose(A),aux(:,count)-q(2:4,count2))
        write(*,'(13F8.3)') x_p(:,count,count2)
      end do
    end do
  end subroutine x_position_body

  subroutine initial_condit (initial_conditions, NumberBodies, q, dq)
    implicit none
    real*8 , allocatable, intent(out):: q(:,:), dq(:,:)
    real*8 , allocatable, intent(in):: initial_conditions(:,:)
    integer , intent(in):: NumberBodies
    integer :: count

    allocate(q(8,NumberBodies+1))
    allocate(dq(8,NumberBodies+1))

    ! Setting initial conditions initial_conditions(13,NumberBodies)
    q = 0d0
    dq = 0d0
    q(1:8,:) = initial_conditions(1:8,:)       !Position
    dq(1,:) = initial_conditions(1,:)          !Velocities
    dq(2:4,:) = initial_conditions(8:10,:) 
    
    do count = 1, NumberBodies + 1
      q(5:8,count) = EulerParam(initial_conditions(5:7,count) )       !Euler parameters
      dq(5:8,count) = dEulerParam(initial_conditions(5:7,count) ,initial_conditions(11:13,count) )      !Euler parameter velocities
    end do
  end subroutine

  subroutine read_initial (iter_max, delta_t, NumberBodies, size_M, constrains, initial_conditions, &
  fixed_coordinates, inertia_tensors, mass)
  
    implicit none

    integer  :: ierr, iter_read, max_iter_read = 1000, count

    integer, intent(out) ::  iter_max
    real*8, intent(out)    :: delta_t

    integer, intent(out)  :: NumberBodies, size_M

    integer, allocatable, intent(out)  :: constrains(:,:)

    real*8 , allocatable, intent(out)  :: initial_conditions(:,:)

    real*8 , allocatable, intent(out) :: fixed_coordinates(:,:)

    real*8 , allocatable, intent(out) :: inertia_tensors(:,:)

    real*8 , allocatable, intent(out) :: mass(:,:)

    open (unit=0,file="initial_input.txt",action="read")

    ! Read parameters solver
    write(*,*) "____________________Numerical parameters___________________"
    do iter_read = 0, max_iter_read
      read (0,*,iostat=ierr) delta_t, iter_max
      if (ierr == 0) then
        write(*,*) delta_t, iter_max
        exit
      end if
    end do

    ! Read number of bodies
    write(*,*) "______________________Number of bodies______________________"
    do iter_read = 0, max_iter_read
      read (0,*,iostat=ierr) NumberBodies
      if (ierr == 0) then
        write(*,*) NumberBodies
        exit
      end if
    end do

    allocate(constrains(8,NumberBodies))
    allocate(initial_conditions(13,NumberBodies+1))
    allocate(fixed_coordinates(8,NumberBodies+1))
    allocate(inertia_tensors(4,NumberBodies*3))
    allocate(mass(2,NumberBodies))

    ! Read constrains
    write(*,*) "__________________________Constrains__________________________"
    count = 1
    do iter_read = 0, max_iter_read
      read (0,*,iostat=ierr) constrains(:,count)
      if (ierr == 0) then
        write(*,*) constrains(:,count)
        count = count + 1
      end if
      if (count == NumberBodies + 1) then
        exit
      end if
    end do

    size_M = NumberBodies*8 + sum(constrains(3:,:))
    write(*,*) "Size M matrix", size_M, "x", size_M

    ! Read initial conditions
    write(*,*) "______________________Initial conditions______________________"
    count = 2
    initial_conditions = 0
    do iter_read = 0, max_iter_read
      read (0,*,iostat=ierr) initial_conditions(:,count)
      if (ierr == 0) then
        write(*,'(13F8.3)') initial_conditions(:,count)
        count = count + 1
      end if
      if (count == NumberBodies + 2) then
        exit
      end if
    end do


    ! Read Fixed coordinates slave bodies
    write(*,*) "______________________Fixed coordinates________________________"
    count = 2
    fixed_coordinates = 0
    if (sum(constrains(3:,:)) /= 0 ) then
      do iter_read = 0, max_iter_read
        read (0,*,iostat=ierr) fixed_coordinates(:,count)
        if (ierr == 0) then
          write(*,'(13F8.3)') fixed_coordinates(:,count)
          count = count + 1
        end if
        if (count == NumberBodies + 2) then
          exit
        end if
      end do
    end if

    ! Read Inertia tensors
    write(*,*) "___________________Inertia tensors____________________"
    count = 1
    do iter_read = 0, max_iter_read
      read (0,*,iostat=ierr) inertia_tensors(:,count)
      if (ierr == 0) then
        write(*,'(13F8.3)') inertia_tensors(:,count)
        count = count + 1
      end if
      if (count == 3*NumberBodies+1) then
        exit
      end if
    end do


    ! Read Mass
    write(*,*) "___________________Mass____________________"
    count = 1
    do iter_read = 0, max_iter_read
      read (0,*,iostat=ierr) mass(:,count)
      if (ierr == 0) then
        write(*,'(13F8.3)') mass(:,count)
        count = count + 1
      end if
      if (count == NumberBodies+1) then
        exit
      end if
    end do

    close(0)

  end subroutine read_initial


  subroutine read_iter(Q_g, ddconstrains_new, NumberBodies,q)
    implicit none
    integer  :: ierr, iter_read, max_iter_read = 1000, count
    real*8, allocatable, intent(out) :: Q_g(:,:), ddconstrains_new(:,:)
    real*8, allocatable :: F_inert(:,:)
    integer, intent(in) :: NumberBodies
    real*8, allocatable, intent(in) :: q(:,:)

    open (unit=0,file="iter_input.txt",action="read")

    allocate(Q_g(8,NumberBodies))
    allocate(F_inert(7,NumberBodies))
    allocate(ddconstrains_new(8,NumberBodies))

    ! Generalized Forces
    !write(*,*) "_________________Generalized forces__________________"
    count = 1
    do iter_read = 0, max_iter_read
      read (0,*,iostat=ierr) F_inert(:,count)
      if (ierr == 0) then
        call inertial_general_F (Q_g(:,count),F_inert(:,count),q(5:8,count+1))
        !write(*,'(13F8.3)') Q_g(:,count)
        count = count + 1
      end if
      if (count == NumberBodies+1) then
        exit
      end if
    end do

    ! Constrains controller
    !write(*,*) "_______________Constrains controller________________"
    count = 1
    do iter_read = 0, max_iter_read
      read (0,*,iostat=ierr) ddconstrains_new(:,count)
      if (ierr == 0) then
        !write(*,'(13F8.3)') ddconstrains_new(:,count)
        count = count + 1
      end if
      if (count == NumberBodies+1) then
        exit
      end if
    end do

    close(0)

  end subroutine read_iter

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
  
  function fM(e, mRR, m00)
    implicit none
      real*8, intent (in) :: e(4), mRR(3,3), m00(3,3)
      real*8 :: fM(8,8), G_(3,4)

      fM = 0d0
      G_ = fG_ (e(1),e(2),e(3),e(4))
      fM(1:3,1:3) = mRR
      fM(4:7,4:7) = MATMUL(MATMUL(TRANSPOSE(G_),m00),G_)
      fM(8,4) = e(1)
      fM(8,5) = e(2)
      fM(8,6) = e(3)
      fM(8,7) = e(4)
      fM(4,8) = e(1)
      fM(5,8) = e(2)
      fM(6,8) = e(3)
      fM(7,8) = e(4)

  end function fM

  function fMconst(x,e,indices)
    implicit none
      real*8, intent (in) :: x(3), e(4)
      integer, intent (in) :: indices(:)
      real*8 :: u(3,3), A(3,3), G(3,4), G_(3,4), auxfMconst(6,7) 
      real*8, allocatable :: fMconst(:,:)
      integer :: dim, index, count

      auxfMconst = 0d0

      G_ = fG_ (e(1),e(2),e(3),e(4))
      G = fG (e(1),e(2),e(3),e(4))
      A = fA (e(1),e(2),e(3),e(4))
      u = fu (x(1),x(2),x(3))
      
      dim = 0
      dim = sum(indices)
      allocate(fMconst(dim,7))
      fMconst = 0d0

      auxfMconst(1,1) = 1d0
      auxfMconst(2,2) = 1d0
      auxfMconst(3,3) = 1d0
      auxfMconst(1:3,4:7) = -MATMUL(MATMUL(A,u),G_)
      auxfMconst(4:6, 4:7) = G

      count = 1
      do index = 1,6
        if (indices(index) == 1) then
          fMconst(count,:) = auxfMconst(index,:)
          count = count + 1
        end if
      end do

  end function fMconst


  function fb(x,e,de,Q,ddc,m00)
    implicit none
      real*8, intent (in) :: x(3), e(4), de(4)
      real*8 :: ddc
      real*8 :: fb(8), u(3,3), G_(3,4), Q(7), dG_(3,4), m00(3,3)

      fb = 0
      G_ = fG_ (e(1),e(2),e(3),e(4))
      dG_ = dfG_ (de(1),de(2),de(3),de(4))
      u = fu (x(1),x(2),x(3))

      fb(1:7) = Q
      fb(4:7) = fb(4:7) - 2d0*MATMUL(MATMUL(MATMUL(TRANSPOSE(dG_),m00),G_), de)
      fb(8) = -DOT_PRODUCT(de,de) + 0.5*ddc
  end function fb


  function fbconst(x,e,de,indices)
    implicit none
    real*8, intent (in) :: x(3), e(4), de(4)
    real*8 ::  auxfbconst(6), u(3,3), dA(3,3), G_(3,4), dG_(3,4)
    real*8, allocatable ::  fbconst(:)
    integer, intent (in) :: indices(:)
    integer :: dim, index, count

    auxfbconst = 0
    G_ = fG_ (e(1),e(2),e(3),e(4))
    dA = dfA (e(1),e(2),e(3),e(4),de(1),de(2),de(3),de(4))
    u = fu (x(1),x(2),x(3))

    dim = 0
    dim = sum(indices)
    allocate(fbconst(dim))
    fbconst = 0d0
    auxfbconst(1:3) = MATMUL(MATMUL(MATMUL(dA,u),G_), de)
    auxfbconst(4:6) = 0 

    count = 1
    do index = 1,6
      if (indices(index) == 1) then
        fbconst(count) = auxfbconst(index)
        count = count + 1
      end if
    end do
  end function fbconst

  function fbcontroller(indices, a_pres)
    implicit none
    real*8, intent (in) :: a_pres(6)
    real*8 ::  auxfbcontroller(6)
    real*8, allocatable ::  fbcontroller(:)
    integer, intent (in) :: indices(:)
    integer :: dim, index, count

    auxfbcontroller = a_pres

    dim = 0
    dim = sum(indices)
    allocate(fbcontroller(dim))
    fbcontroller= 0d0

    count = 1
    do index = 1,6
      if (indices(index) == 1) then
        fbcontroller(count) = auxfbcontroller(index)
        count = count + 1
      end if
    end do
  end function fbcontroller

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
