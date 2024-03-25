module common
    implicit none
    integer, public, parameter :: double = selected_real_kind(13)
PUBLIC :: initial, wolff, chose_random_seed_spin, init_random_seed, look_at_neigh
PUBLIC ::  flip_spin, coor
real (kind = double) , public :: b,E,m, padd, j_int=1.
INTEGER , public , dimension ( :, :, : ), ALLOCATABLE :: spin, cluster, neigh
integer, public :: num,L,nmcs,nequil, seed_spin, cx, cy, cz, n_clu, c, c2, imcs
INTEGER , public , dimension ( : ), ALLOCATABLE :: coo
contains


subroutine initial()

integer ::  x, y, z, up, right, front, sums
real (kind =double) :: rnd
     !print *, "linear dimension of lattice L ="
    !read *, L
!L = 8
!b = 2
j_int = 1
    allocate(spin(L,L,L))
    allocate(cluster(L,L,L))
    allocate(neigh(L,L,L))
    
    !print *, "reduced temperature T ="
    !read *, b
    num = L*L*L
    allocate(coo(3*num))
    padd = 1 - exp(-2*j_int/b)
    !print *, padd
    !print *, "# MC steps per spin for equilibrium ="
    !read *, nequil
    !print *, "# MC steps per spin for averages ="
    !read *, nmcs
    call init_random_seed()
    n_clu = 0
    
    !  random initial configuration
    !  compute initial magnetization
    do z =1,L
    do y = 1,L
       do x = 1,L
          call random_number(rnd)
          if (rnd < 0.5) then
             spin(x,y, z) = 1
          else
             spin(x,y,z) = -1
          end if
          
       end do
    end do
end do
    !  compute initial energy
    E = 0.0_double
    do z=1,L
        if (z == L) then
            front = 1
            
         else
            front = z + 1
         end if
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
          sums = spin(x,up, z) + spin(right,y,z) + spin(x,y,front)
! calculate the initial energy summing all over pairs
! (gor a given spin, consider only the up NN and the right NN
! - NOT the down and the left NN - : each interaction is counted once
          E = E - spin(x,y,z)*sums
       end do
    end do
end do
    end subroutine initial




subroutine wolff()

integer :: i, ncx, ncy, ncz, io
integer ::  x, y, z, up, right, front, sums
i = 0
  call chose_random_seed_spin()
    do 
       
  call look_at_neigh() 
  io = 0
        call coor(1+3*i, ncx, ncy, ncz, io) 
        cx = ncx
        cy = ncy
        cz = ncz
        if ( i == n_clu -1) exit
        i = i+1
    end do
    call cluster_size()
    call flip_spin()
    m = 0
    do z=1,L
        do y =1,L
            do x =1,L
                m = m + spin(x,y,z) 
            end do
        end do
    end do

    E = 0.0_double
    do z=1,L
        if (z == L) then
            front = 1
            
         else
            front = z + 1
         end if
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
          sums = spin(x,up, z) + spin(right,y,z) + spin(x,y,front)
          E = E - spin(x,y,z)*sums
       end do
    end do
end do
  end subroutine wolff


  subroutine data(cum)
    !  accumulate data after every Monte Carlo step per spin
    real (kind = double), dimension(5), intent (inout) :: cum
    real (kind = double) :: eave,e2ave,mave,m2ave, esigma2, msigma2
    cum(1) = cum(1) + E
    cum(2) = cum(2) + E*E
    cum(3) = cum(3) + m
    cum(4) = cum(4) + m*m
    cum(5) = cum(5) + abs(m)
    eave     = cum(1)/real(num)/real(imcs)    ! to avoid interger overflow
    e2ave    = cum(2)/real(num*num)/real(imcs)
    mave     = cum(3)/real(num)/real(imcs)
    m2ave    = cum(4)/real(num*num)/real(imcs)
    esigma2 = e2ave - (eave*eave)
    msigma2 = m2ave - (mave*mave)
    !write(51, *) imcs,  mave,  msigma2/(real(b)*real(num))
  end subroutine data

  subroutine output(cum)
    real (kind = double), dimension(5), intent (inout) :: cum
    real (kind = double) :: eave,e2ave,mave,m2ave,abs_mave, esigma2, msigma2
    eave     = cum(1)/real(num)/real(nmcs)    ! to avoid interger overflow
    e2ave    = cum(2)/real(num*num)/real(nmcs)
    mave     = cum(3)/real(num)/real(nmcs)
    m2ave    = cum(4)/real(num*num)/real(nmcs)
    abs_mave = cum(5)/real(num)/real(nmcs)
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
    write(31, *) b, esigma2/(real(b)*real(b)*real(num))
   ! df(1)=df(2)
    
  end subroutine output

subroutine chose_random_seed_spin()


    real ( kind = double ) :: rx, ry, rz 
    cluster = 0
    neigh = 0
    coo = 0
    call random_number( rx)
    call random_number( ry)
    call random_number( rz)
    cx=nint(rx*(L-1)) +1
    cy=nint(ry*(L-1)) +1
    cz=nint(rz*(L-1)) +1 
    seed_spin = spin(cx, cy,cz)   
    !print *, 'seed', cx, cy,cz, seed_spin
    cluster(cx, cy,cz) = 1
    n_clu = 1
    c = 0
    call coor(n_clu, cx, cy, cz, 1)

end subroutine chose_random_seed_spin

subroutine look_at_neigh()

    integer :: i, nx, ny,nz




do i=0,2
    nx = cx -1+i
    if (nx == L+1)  nx = 1
    if (nx == 0)  nx = L
    if ( cluster(nx, cy,cz) == 0 ) call add(nx, cy,cz)
    !print *, 'nx', nx, cy,cz, spin(nx, cy,cz)
end do

do i=0,2
    nz = cz -1+i
    if (nz == L+1)  nz = 1
    if (nz == 0)  nz = L
    
    if (cluster(cx, cy,nz) == 0 ) call add(cx, cy,nz)
    !print *, 'nz', cx, cy,nz, spin(cx, cy,nz)
end do

do i=0,2
    ny = cy -1+i
    if (ny == L+1)  ny = 1
    if (ny == 0)  ny = L
    
    if ( cluster(cx, ny,cz) == 0) call add(cx, ny,cz)
    !print *, 'ny' ,cx, ny,cz, spin(cx, ny,cz)
    


end do

    end subroutine look_at_neigh

 subroutine add(nx,  ny, nz)
        
        integer ,intent(in) :: nx,  ny, nz
        real (kind = double) :: rnd
        integer :: ncx, ncy, ncz, io
        
        call random_number(rnd)
        if (spin(nx,ny,nz) == seed_spin) then
            !print *, 'add', nx,ny,nz, spin(nx,ny,nz)
            !print *, rnd, padd
             if (rnd <= padd) then
                !print *, rnd, padd
                 cluster(nx,ny,nz) = 1
                 n_clu=n_clu+1
                 !print *, 'n_clu' ,n_clu
                 ncx = nx
                 ncy = ny
                 ncz = nz
                 io = 1
                 call coor(1+ 3*(n_clu-1), ncx,ncy,ncz, io)
                 
            end if
           
        end if 



        end subroutine add



subroutine coor(i, nx,  ny, nz, io)

    INTEGER, INTENT(IN) :: i 
INTEGER,intent(inout) ::  nx,  ny, nz
integer, INTENT(IN) :: io

if ( io == 1 ) then
    coo(i) = nx
    coo(i+1) = ny
    coo(i+2) = nz
else
     nx = coo(i) 
     ny = coo(i+1) 
     nz = coo(i+2) 
end if


end subroutine coor






subroutine flip_spin()

    INTEGER :: i, nxc, ncy, ncz, io
    io = 0
do i=0, n_clu-1
    call coor(1 +3*i, nxc, ncy, ncz, io) !!!!!!!!!!!
    spin(nxc, ncy, ncz) = spin(nxc, ncy, ncz)* (-1)
end do

end subroutine flip_spin

subroutine cluster_size()

    integer :: i, k, v
    do i =1,L
        do k = 1,L
            do v = 1,L
        if ( cluster(i, k, v) > 0 ) then
            c = c+1
        end if
    end do
end do
    end do


end subroutine cluster_size











subroutine init_random_seed()
    use iso_fortran_env, only: int64
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid
    integer(int64) :: t
  
    call random_seed(size = n)
    allocate(seed(n))
    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(un) seed
       close(un)
    else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       call system_clock(t)
       if (t == 0) then
          call date_and_time(values=dt)
          t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_int64 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
       end if
       pid = getpid()
       t = ieor(t, int(pid, kind(t)))
       do i = 1, n
          seed(i) = lcg(t)
       end do
    end if
    call random_seed(put=seed)
  contains
    ! This simple PRNG might not be good enough for real work, but is
    ! sufficient for seeding a better PRNG.
    function lcg(s)
      integer :: lcg
      integer(int64) :: s
      if (s == 0) then
         s = 104729
      else
         s = mod(s, 4294967296_int64)
      end if
      s = mod(s * 279470273_int64, 4294967291_int64)
      lcg = int(mod(s, int(huge(0), int64)), kind(0))
    end function lcg
  end subroutine init_random_seed


end module common




program isingwolff
    
    use common

    integer :: p
    real (kind = double), dimension(5) :: cum
    nequil = 300
    nmcs = 10000
    do L = 8,8
    do p = 75, 100
        b = real(p)/20
        !b = 5
   call initial
   cum = 0
   !print *, 'padd', b,padd
   do imcs=1, nequil
    call wolff()
   end do
   m = 0
   do imcs=1, nmcs
   call wolff()

   call data(cum)
   end do

   call output(cum) 
   
    DEALLOCATE(spin)
    DEALLOCATE(cluster)
    DEALLOCATE(neigh)
    DEALLOCATE(coo)
    end do
end do
  end program isingwolff
  