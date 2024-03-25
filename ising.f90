!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c ising.f90
!c
!c Metropolis algorithm to calculate <E>, <M>, in the canonical ensemble
!c (fix T,N,V) with a 2D Ising model
!c
!c Here: K_B = 1
!c       J   = 1
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module common

  implicit none
  public :: initial,metropolis,DeltaE
  public :: data,output
  integer, public, parameter :: double = selected_real_kind(13)
  real (kind = double), public :: T,E,M
  integer, public, dimension(:,:,:), allocatable :: spin
  real (kind = double), public, dimension(-20:20) :: w
  integer, public :: N,L,nmcs,nequil
  integer, public :: accept
  real (kind = double), dimension(2), public :: df
contains

  subroutine initial(nequil,cum)
    integer, intent (out) :: nequil
    real (kind = double), dimension(5), intent (out) :: cum
    integer :: x,y,up,right,sums,dE, sizer, z, front
    integer, dimension(:), allocatable :: seed
    real :: rnd
    df=0.0_double
    !L=8
    nequil=300
    nmcs=10000
    !T=2
    !print *, "linear dimension of lattice L ="
    !read *, L
    allocate(spin(L,L,L))
    !print *, "reduced temperature T ="
    !read *, T
    N = L*L*L
    !print *, "# MC steps per spin for equilibrium ="
    !read *, nequil
    !print *, "# MC steps per spin for averages ="
    !read *, nmcs
    call random_seed(sizer)
    allocate(seed(sizer))
    call random_seed(get=seed)
    M = 0.0_double
    !  random initial configuration
    !  compute initial magnetization
    do y = 1,L
       do x = 1,L
         do z = 1,L
          call random_number(rnd)
          if (rnd < 0.5) then
             spin(x,y, z) = 1
          else
             spin(x,y, z) = -1
          end if
          M = M + spin(x,y,z)
         end do
       end do
    end do
    !  compute initial energy
    E = 0.0_double
    do y = 1,L
       !  periodic boundary conditions
       if (y == L) then
          up = 1
          
       else
          up = y + 1
       end if
       do x = 1,L
          if (x == L) then
             right = 1
             
          else
             right = x + 1
          end if
          do z = 1,L
            if (x == L) then
               front = 1
               
            else
               front = x + 1
            end if
          sums = spin(x,up, z) + spin(right,y, z) + spin(x,y, front)
! calculate the initial energy summing all over pairs
! (gor a given spin, consider only the up NN and the right NN
! - NOT the down and the left NN - : each interaction is counted once
          E = E - spin(x,y, z)*sums
         end do
       end do
    end do
!
! calculate the transition probability according 
! to the Boltzmann distribution (exp(-deltaE/KT).
! Choosing the interaction parameter  J=1, ***IN CASE OF P.B.C.***
! possible energy variations per spin flip are -8,-4,0,+4,+8: 
    do dE = -8,8,4
       w(dE) = exp(-dE/T)
    end do
    accept = 0
    cum = 0.0_double
  end subroutine initial

  subroutine metropolis()
    !  one Monte Carlo step per spin
    integer :: ispin,x,y,dE, z
    real :: rnd
    do ispin = 1,N
       !     random x and y coordinates for trial spin
       call random_number(rnd)
       x = int(L*rnd) + 1
       call random_number(rnd)
       y = int(L*rnd) + 1
       call random_number(rnd)
       z = int(L*rnd) + 1
       dE = DeltaE(x,y, z)
       call random_number(rnd)
       if (rnd <= w(dE)) then
          spin(x,y, z) = -spin(x,y, z)
          accept = accept + 1
          M = M + 2*spin(x,y, z)  ! factor 2 is to account for the variation:
          E = E + dE           ! (-(-)+(+))
       end if
    end do
  end subroutine metropolis

  function DeltaE(x,y,z) result (DeltaE_result)
    !  periodic boundary conditions
    integer, intent (in) :: x,y,z
    integer :: DeltaE_result
    integer :: left
    integer :: right
    integer :: up
    integer :: down, front, beh
    if (x == 1) then
       left = spin(L,y, z)
       right = spin(2,y, z)
       !left = 0
    else if (x == L) then
       left = spin(L-1,y, z)
       right = spin(1,y, z)
       !right = 0
    else
       left = spin(x-1,y, z)
       right = spin(x+1,y, z)
    end if
    if (y == 1) then
       up = spin(x,2, z)
       down = spin(x,L, z)
       !down = 0
    else if (y == L) then
       up = spin(x,1, z)
       down = spin(x,L-1, z)
       !up = 0 
    else
       up = spin(x,y+1, z)
       down = spin(x,y-1, z)
    end if
    if (z == 1) then
      front = spin(x,y, 2)
      beh = spin(x,y, L)
      !beh = 0
   else if (z == L) then
      front = spin(x,y, 1)
      beh = spin(x,y, L-1)
      !front = 0 
   else
      front = spin(x,y, z+1)
      beh = spin(x,y, z-1)
   end if
    DeltaE_result = 2*spin(x,y,z)*(left + right + up + down + front + beh)
! also here the factor 2 is to account for the variation
  end function DeltaE

  subroutine data(cum, imcs)
    !  accumulate data after every Monte Carlo step per spin
    real (kind = double), dimension(5), intent (inout) :: cum
    integer, intent(in) :: imcs
    real (kind = double) :: eave,e2ave,mave,m2ave, esigma2, msigma2
    cum(1) = cum(1) + E
    cum(2) = cum(2) + E*E
    cum(3) = cum(3) + M
    cum(4) = cum(4) + M*M
    cum(5) = cum(5) + abs(M)
    eave     = cum(1)/real(N)/real(imcs)    ! to avoid interger overflow
    e2ave    = cum(2)/real(N*N)/real(imcs)
    mave     = cum(3)/real(N)/real(imcs)
    m2ave    = cum(4)/real(N*N)/real(imcs)
    esigma2 = e2ave - (eave*eave)
    msigma2 = m2ave - (mave*mave)
    !write(11, *) imcs, eave, mave, esigma2/(real(T)*real(T)*real(N)), msigma2/(real(T)*real(N))
    !write(51, *) imcs, mave,  msigma2/(real(T)*real(N))
  end subroutine data

  subroutine output(cum)
    real (kind = double), dimension(5), intent (inout) :: cum
    real (kind = double) :: eave,e2ave,mave,m2ave,abs_mave, esigma2, msigma2
    real :: acceptance_prob
    acceptance_prob = real(accept)/real(N)/real(nmcs+nequil)
    eave     = cum(1)/real(N)/real(nmcs)    ! to avoid interger overflow
    e2ave    = cum(2)/real(N*N)/real(nmcs)
    mave     = cum(3)/real(N)/real(nmcs)
    m2ave    = cum(4)/real(N*N)/real(nmcs)
    abs_mave = cum(5)/real(N)/real(nmcs)
    esigma2 = e2ave - (eave*eave)
    msigma2 = m2ave - (abs_mave*abs_mave)
   ! print *, "temperature                =", T
   ! print *, "acceptance probability     =", acceptance_prob
   ! print *, "mean energy per spin                =", eave
    !print *, "mean squared energy per spin        =", e2ave
   ! print *, "mean magnetization  per spin        =", mave
    !print *, "mean squared magnetization per spin =", m2ave
    !print *, "mean |magnetization| per spin       =", abs_mave
   ! print *, "heat capacity per spin     =", esigma2/(real(T)*real(T)*real(N))
   ! print *, "magnetic susceptibility per spin     =", msigma2/(real(T)*real(N))
   ! df(2)=esigma2/(real(T)*real(T)*real(N))
    write(10, *) T, esigma2/(real(T)*real(T)*real(N))
   ! df(1)=df(2)
    
  end subroutine output

end module common

program ising
  !  metropolis algorithm for the ising model on a square lattice
  use common
  integer :: imcs,ispin,jspin, i
  real (kind = double), dimension(5) :: cum
  do L = 8,8
  do i=75, 100
    T=real(i)/20
   !T = 5
     call initial(nequil,cum)
     !  equilibrate system
     do imcs = 1,nequil
        call metropolis()
        !do jspin = 1,L
         !     do ispin = 1,L
         !        if(spin(ispin,jspin)==1)write(20,*)ispin,jspin
         !       
         !        if(spin(ispin,jspin)==-1)write(40,*)ispin,jspin
          !      
          !    end do
          ! end do
        !  write(20,*)
        !  write(20,*)
        !  write(40,*)
        !  write(40,*)
     end do
     !  accumulate data while updating spins
     
     do imcs = 1,nmcs
        call metropolis()
        call data(cum, imcs)
        !if(mod(imcs,100) == 0) then
          ! do jspin = 1,L
           !   do ispin = 1,L
            !     if(spin(ispin,jspin)==1)write(20,*)ispin,jspin
                
             !    if(spin(ispin,jspin)==-1)write(40,*)ispin,jspin
                
            !  end do
          ! end do
          !write(20,*)
          !write(20,*)
         ! write(40,*)
         ! write(40,*)
        !end if
     end do
  
     call output(cum)
     deallocate(spin)
  end do
end do
! write the coordinates of spins up and down on files for plotting
  !open(unit=8,file='ising-up.dat',status='replace')
  !open(unit=9,file='ising-down.dat',status='replace')
  !do jspin = 1,L
   !  do ispin = 1,L
    !   if(spin(ispin,jspin)==1)write(8,*)ispin,jspin
     !   if(spin(ispin,jspin)==-1)write(9,*)ispin,jspin
     !end do
 ! end do
  !close(8)
  !close(9)
  
  !deallocate(spin)
end program ising
